""" Helper functions """

__author__ = "Vince Reuter"
__email__ = "vince.reuter@gmail.com"


def get_epi_cmd(jar, readsfile, sitesfile, outdir, min_rlen, min_qual,
    strand_method, rrbs_fill, mem_gig, context="CG", cores=1):
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
    outdir : str
        Path to folder for output files.
    min_rlen : int
        Minimum number of bases aligned for a read to be used in analysis.
    min_qual : int
        Minimum base call quality at a site and neighbor site(s) for it to
        be used in analysis.
    strand_method : str
        Name of strategy to determine read orientation; 'tag' or 'flag'
    rrbs_fill : int
        Number of bases at read end to ignore due to RRBS "fill-in"
    mem_gig : int
        Number of gigabytes of memory.
    context : str
        Methylation context (sense strand, e.g. 'CG' for typical mammalian analysis)
    cores : int
        Number of cores to use for processing

    Returns
    -------
    str
        Base of command for epiallele processing.

    """

    import os

    contexts = ["C", "CG"]

    def mem_err_msg():
        return "Mem spec must be nonnegative integer; got {}".format(mem_gig)

    problems = []

    pos_int_vals = {"Memory": mem_gig, "Length": min_rlen, "Quality": min_qual}
    for label, value in pos_int_vals.items():
        valid = False
        try:
            valid = int(value) >= 0
        except (TypeError, ValueError):
            pass
        if not valid:
            problems.append("Did not get nonnegative value -- {} = {}".format(label, value))

    try:
        mem_gig = int(mem_gig)
    except (TypeError, ValueError):
        problems.append(mem_err_msg())
    else:
        if mem_gig < 0:
            problems.append(mem_err_msg())

    problems.extend(
        ["{} isn't a file: {}".format(ft, fp)
         for ft, fp in zip(["JAR", "Reads", "Sites"], [jar, readsfile, sitesfile])
         if not os.path.isfile(fp)])

    if context not in contexts:
        problems.append("Invalid context ({}); choose one: {}".format(context, ", ".join(contexts)))

    try:
        if int(cores) < 1:
            problems.append("Too few cores: {}".format(cores))
    except:
        problems.append("Invalid cores count: {}".format(cores))

    if problems:
        raise Exception("Problems: {}".format(", ".join(problems)))
    return "java -Xmx{m}g -jar {j} --minBaseQuality {q} --minReadLength {rl} --context {ctx} --rrbsFill {base_fill} --cores {cores} --strandMethod {sm} -O {o} {r} {s}".format(
        m=mem_gig, j=jar, rl=min_rlen, q=min_qual, ctx=context, base_fill=rrbs_fill, sm=strand_method, o=outdir, r=readsfile, s=sitesfile, cores=cores)
