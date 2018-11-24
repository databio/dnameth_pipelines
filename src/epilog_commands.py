""" Helper functions and data types """

import os
from helpers import EpilogTarget

__author__ = "Vince Reuter"
__email__ = "vince.reuter@gmail.com"

__all__ = [
    "get_epilog_full_command", "get_epilog_union_command",
    "get_epilog_union_calls_command", "get_epilog_union_epis_command",
    "get_epilog_strand_merge_command", "get_epilog_epistats_command",
    "get_proc_stat_log", "make_main_epi_cmd", "run_main_epi_pipe"]


# Downstream processing constants
SINGLE_SITES_DATA_TYPE_NAME = "sites"
EPIALLELES_DATA_TYPE_NAME = "epialleles"
SUFFIX_BY_TYPE = {SINGLE_SITES_DATA_TYPE_NAME: "calls", EPIALLELES_DATA_TYPE_NAME: "epialleles"}


def make_main_epi_cmd(
    epiconf, prog_spec, readsfile, sitesfile, outdir,
    rrbs_fill, context, epis=True, process_logfile=None):
    """
    Simplified version of the builder of the main epilog processing command.

    :param pypiper.AttributeDict epiconf: Methylation downstream analysis
        configuration/parameterization
    :param helpers.ProgSpec prog_spec: JAR, memory allocation, and cores count
    :param str readsfile: Sorted, aligned BAM
    :param str sitesfile: Gzipped tabix-indexed collection of genome positions
        at which to look for methylation
    :param str outdir: Path to folder in which to place output
    :param int rrbs_fill: Number of bases to ignore for "RRBS fill-in"
    :param str context: Sequence context
    :param bool epis: Whether to examine epialleles
    :param str process_logfile: Path to file to which to write epilog
        basic processing performance numbers
    :return str, Sequence of str: Command to run main epilog processing, and
        sequence of targets (files) that it should produce
    """
    return get_epilog_full_command(prog_spec, readsfile, sitesfile, outdir,
        min_rlen=epiconf.read_length_threshold, min_qual=epiconf.qual_threshold,
        strand_method=epiconf.strand_method, rrbs_fill=rrbs_fill,
        context=context, epis=epis, halt="union", process_logfile=process_logfile)


def get_epilog_full_command(prog_spec, readsfile, sitesfile, outdir,
    min_rlen, min_qual, strand_method, rrbs_fill, context, epis=True,
    halt=None, strand_specific=False, no_epi_stats=False, process_logfile=None):
    """
    Create base for epiallele processing command.

    Parameters
    ----------
    prog_spec : ProgSpec
        Bundle of JAR, memory allocation, and number of cores
    readsfile : str
        Path to sorted, aligned BAM with reads to analyze.
    sitesfile : str
        Path to gzipped, tabix-indexed file with sites to analyze.
    outdir : str
        Path to output folder for single-site (and epiallele, as desired) output.
    min_rlen : int
        Minimum number of bases aligned for a read to be used in analysis.
    min_qual : int
        Minimum base call quality at a site and neighbor site(s) for it to
        be used in analysis.
    strand_method : str
        Name of strategy to determine read orientation; 'tag' or 'flag'
    rrbs_fill : int
        Number of bases at read end to ignore due to RRBS "fill-in"
    context : str
        Methylation context (sense strand, e.g. 'CG' for typical mammalian analysis)
    epis : bool
        Produce epiallele results in addition to the more single-site-oriented output
    halt : str, optional
        Name of processing stage after which to halt; if omitted, run as much as possible
    strand_specific : bool, default False
        Indicate no strand merger is desired.
    no_epi_stats : bool, default True
        Skip epiallele diversity/heterogeneity statistics
    process_logfile : str, optional
        Path to file for epiallele processing performance statistics

    Returns
    -------
    str, str or list of str
        Command for main epilog processing, and a pypiper "target"
        (path to calls file, or that and epiallele path if applicable)

    """

    import os

    contexts = ["C", "CG"]

    problems = []

    pos_int_vals = {"Length": min_rlen, "Quality": min_qual}
    for label, value in pos_int_vals.items():
        valid = False
        try:
            valid = int(value) >= 0
        except (TypeError, ValueError):
            pass
        if not valid:
            problems.append("Did not get nonnegative value -- {} = {}".format(label, value))

    problems.extend(["Missing {} file: {}".format(ft, fp) for ft, fp in
        zip(["reads", "sites"], [readsfile, sitesfile]) if not os.path.isfile(fp)])

    if context not in contexts:
        problems.append("Invalid context ({}); choose one: {}".format(context, ", ".join(contexts)))

    if problems:
        raise Exception("Problems: {}".format(", ".join(problems)))

    name_ss_file = "all_calls.txt"

    def get_outpath(fn):
        return os.path.join(outdir, fn)

    single_sites_file = get_outpath(name_ss_file)

    cmd = "{b} --minBaseQuality {q} --minReadLength {rl} --context {ctx} --rrbsFill {base_fill} --cores {cores} --strandMethod {sm} -O {o} {r} {s}".format(
        b=prog_spec.get_command_base(), rl=min_rlen, q=min_qual, ctx=context,
        base_fill=rrbs_fill, sm=strand_method, o=single_sites_file, r=readsfile, s=sitesfile, cores=prog_spec.cores)

    if epis:
        epis_file = get_outpath("all_epialleles.txt")
        cmd += " --outputEpialleles {}".format(epis_file)
    else:
        epis_file = None
    if process_logfile:
        cmd += " --processLogfile {}".format(process_logfile)
    if halt:
        cmd += " --through {}".format(halt)
    if strand_specific:
        cmd += " --strandSpecific"
    if no_epi_stats:
        cmd += " --noEpiStats"

    # TODO: though this is not going to be the encouraged route, while/if it's to be provided, consider the downstream file(s) as targets.
    # TODO: beware, though, of the effect on the "main-only" function that calls into this. Its targets are the main files.
    return cmd, EpilogTarget(single_sites_file=single_sites_file, epis_file=epis_file)


def get_epilog_union_command(prog_spec, data_type, folder, output=None):
    """
    Create command with which to union per-chromosome (or other block/chunk) results files.

    :param ProgSpec prog_spec: Bundle of JAR, memory, and cores count
    :param str folder: Path to the folder with files to union
    :param str data_type: Name indicating type of input data; choose 'sites' or 'epialleles'
    :param str output: Path to output file; optional, can be derived from data type
    :return str, str: Command and output path
    """
    data_type = _validate_data_type(data_type)
    output = output or os.path.join(folder, "all_{}.txt".format(SUFFIX_BY_TYPE[data_type]))
    base = prog_spec.get_command_base(("epiallele", "UnionCallsFiles"))
    cmd = "{b} {dt} -I {i} -O {o}".format(b=base, dt=data_type, i=folder, o=output)
    return cmd, output


def get_epilog_union_calls_command(prog_spec, folder, output=None):
    return get_epilog_union_command(prog_spec, data_type=SINGLE_SITES_DATA_TYPE_NAME, folder=folder, output=output)


def get_epilog_union_epis_command(prog_spec, folder, output=None):
    return get_epilog_union_command(prog_spec, data_type=EPIALLELES_DATA_TYPE_NAME, folder=folder, output=output)


def get_epilog_strand_merge_command(
    prog_spec, data_file, data_type, outfolder=None, gff=False, out_ext=".merged.txt"):
    """
    Create command to merge data for the same sites but from different strand

    :param ProgSpec prog_spec: Bundle of JAR, memory, and cores count
    :param str data_file: Path to file with data to merge
    :param str data_type: Name for type of data being merged
    :param str outfolder: Path to output folder; if unspecified, same as input
    :param bool gff: Whether the data to merge are stored as GFF
    :param str out_ext: Extension for resulting file
    :return str, str: Command and path to output file
    """
    dt = _validate_data_type(data_type)
    infolder, infile = os.path.split(data_file)
    outfolder = outfolder or infolder
    outname = os.path.splitext(infile)[0] + out_ext
    output = os.path.join(outfolder, outname)
    base = prog_spec.get_command_base(("analysis", "CombineStrands"))
    cmd = "{b} -I {i} -O {o} -F {f} -T {t}".format(b=base, i=data_file, o=output, f="gff" if gff else "delim", t=dt)
    return cmd, output


def get_epilog_epistats_command(prog_spec, infile, outfile, stranded=None, gff=False):
    """
    Create command to calculate epiallele heterogeneity statistics.

    :param ProgSpec prog_spec: Bundle of JAR, memory, and cores count
    :param str infile: Path to input (epiallele calls) file
    :param str outfile: Path to file to which to write output.
    :param bool stranded: Indicate that input records have strand column
    :param bool gff: Indicate that input data is stored as GFF
    :return str, str: Command and path to output file
    """
    base = prog_spec.get_command_base(("analysis", "AnalyzeEpihet"))
    cmd = "{b} -I {i} -O {o}".format(b=base, i=infile, o=outfile)
    if gff:
        cmd += " --gff"
        if stranded is False:
            print("WARNING: lack of strand is illegal for GFF; ignoring no-strand setting")
    elif stranded:
        cmd += " --hasStrandCol"
    return cmd, outfile


def get_proc_stat_log(folder):
    """
    From epilog output folder, get path to file for processing performance.

    :param str folder: Path to folder in which to place the processing stats file
    :return str: Path to file for processing stats
    """
    return os.path.join(folder, "processing_performance.log")


def run_main_epi_pipe(pm, epiconf, prog_spec, readsfile, sitesfile, outdir, rrbs_fill):
    """
    Run downstream methylation analysis for a sample's "real" (non-spike-in) data.

    :param pypiper.PipelineManager pm: Pipeline manager, to run commands
    :param pypiper.AttributeDict epiconf: epilog parameterization
    :param helpers.ProgSpec prog_spec: JAR, memory spec text, and cores count
    :param str readsfile: Path to aligned, dedup, sorted, indexed BAM
    :param str sitesfile: Path to gzipped, tabix-indexed collection of methylation
        positions to consider
    :param str outdir: Path to epilog output folder
    :param int rrbs_fill: Number of bases to disregard on account of RRBS
    """

    from helpers import missing_targets

    epi_main_cmd, epi_main_tgt = make_main_epi_cmd(
        epiconf, prog_spec, readsfile, sitesfile, outdir,
        rrbs_fill=rrbs_fill, context=epiconf.context,
        epis=True, process_logfile=get_proc_stat_log(outdir))
    pm.run(epi_main_cmd, target=epi_main_tgt, lock_name="epilog_main", nofail=True)

    # Proceed with strand merger (if desired) based on the presence of the targets.
    missing = missing_targets(epi_main_tgt)
    if missing:
        print("Missing main epilog target(s): {}".format(", ".join(missing)))
    elif not epiconf.strand_specific:
        pm.timestamp("### Epilog strand merger")
        merge_cmd, merged_epi_tgt = get_epilog_strand_merge_command(
            prog_spec, epi_main_tgt.epis_file, data_type="epialleles")
        pm.run(merge_cmd, merged_epi_tgt, lock_name="epilog_merge_epis", nofail=True)
        merge_cmd, merged_ss_tgt = get_epilog_strand_merge_command(
            prog_spec, epi_main_tgt.single_sites_file, data_type="sites")
        pm.run(merge_cmd, merged_ss_tgt, lock_name="epilog_merge_single", nofail=True)
        epis_file = merged_epi_tgt
    else:
        epis_file = epi_main_tgt.epis_file

    if not epiconf.no_epi_stats:
        def exp_skip(exp):
            print("Due to {}, {}".format(exp, "epiallele statistics cannot be calculated."))
        if missing:
            exp_skip("missing results upstream")
        elif not os.path.isfile(epis_file):
            exp_skip("missing epialleles file ({})".format(epis_file))
        else:
            epilog_stats_target = os.path.join(outdir, "epiallele_statistics.txt")
            epi_stats_cmd, epi_stats_tgt = get_epilog_epistats_command(prog_spec,
                infile=epis_file, outfile=epilog_stats_target, stranded=epiconf.strand_specific)
            pm.run(epi_stats_cmd, epi_stats_tgt, lock_name="epilog_epistats", nofail=True)


def _validate_data_type(dt):
    """ Check that the data type name specification complies with the tool used. """
    try:
         dt = dt.lower()
    except AttributeError:
        raise TypeError("Data type to union must be string, but got {}".format(type(dt)))
    if dt not in SUFFIX_BY_TYPE:
        raise ValueError("Illegal data type name ('{}'); choose from: {}".format(dt, ", ".join(SUFFIX_BY_TYPE.keys())))
    return dt
