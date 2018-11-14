""" Helper functions """

__author__ = "Vince Reuter"
__email__ = "vince.reuter@gmail.com"


def get_epi_cmd(jar, readsfile, sitesfile, outfile, min_rlen, min_qual,
    strand_method, rrbs_fill, memtext,
    context="CG", cores=1, keep_chrom_files=False,
    epis_file=None, process_logfile=None):
    """
    Create base for epiallele processing command.

    Parameters
    ----------
    jar : str
        Path to JAR file for epiallele software.
    readsfile : str
        Path to sorted, aligned BAM with reads to analyze.
    sitesfile : str
        Path to gzipped, tabix-indexed file with sites to analyze.
    outfile : str
        Path to file for site-level methylation call data (main target).
    min_rlen : int
        Minimum number of bases aligned for a read to be used in analysis.
    min_qual : int
        Minimum base call quality at a site and neighbor site(s) for it to
        be used in analysis.
    strand_method : str
        Name of strategy to determine read orientation; 'tag' or 'flag'
    rrbs_fill : int
        Number of bases at read end to ignore due to RRBS "fill-in"
    memtext : str
        Text specification of memory, e.g. 16000m or 4g.
    context : str
        Methylation context (sense strand, e.g. 'CG' for typical mammalian analysis)
    cores : int
        Number of cores to use for processing.
    keep_chrom_files : bool
        Whether the per-chromosome output files from epilog should be retained
        along with the main, merged output files.
    epis_file : str, optional
        Path to file for epiallele observation records; if unspecified, no
        epiallele processing will be performed.
    process_logfile : str, optional
        Path to file for epiallele processing performance statistics

    Returns
    -------
    str
        Base of command for epiallele processing.

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

    problems.extend(
        ["{} isn't a file: {}".format(ft, fp)
         for ft, fp in zip(["JAR", "Reads", "Sites"], [jar, readsfile, sitesfile])
         if not os.path.isfile(fp)])

    if context not in contexts:
        problems.append("Invalid context ({}); choose one: {}".format(context, ", ".join(contexts)))

    try:
        if int(cores) < 1:
            problems.append("Too few cores: {}".format(cores))
    except (TypeError, ValueError):
        problems.append("Invalid cores count: {}".format(cores))

    if problems:
        raise Exception("Problems: {}".format(", ".join(problems)))

    cmd = "java -Xmx{m} -jar {j} --minBaseQuality {q} --minReadLength {rl} --context {ctx} --rrbsFill {base_fill} --cores {cores} --strandMethod {sm} -O {o} {r} {s}".format(
        m=memtext, j=jar, rl=min_rlen, q=min_qual, ctx=context, base_fill=rrbs_fill, sm=strand_method, o=outfile, r=readsfile, s=sitesfile, cores=cores)

    if keep_chrom_files:
        cmd += " --keepChromFiles"
    if epis_file:
        cmd += " --outputEpialleles {}".format(epis_file)
    if process_logfile:
        cmd += " --processLogfile {}".format(process_logfile)

    return cmd
