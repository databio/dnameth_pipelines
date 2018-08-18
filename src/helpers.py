""" Helper functions """

__author__ = "Vince Reuter"
__email__ = "vince.reuter@gmail.com"


def get_epi_cmd(jar, readsfile, sitesfile, outdir, strand_method, rrbs_fill, mem_gig, context="CG"):
    """
    Create base for epiallele processing command.

    Parameters
    ----------
    jar : str
        Path to JAR file for epiallele software.
    mem_gig : int
        Number of gigabytes of memory.

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

    if problems:
        raise Exception("Problems: {}".format(", ".join(problems)))
    return "java -Xmx{}g -jar {} --context {} --rrbsFill {} --strandMethod {} -O {} {} {}".format(
        mem_gig, jar, context, rrbs_fill, strand_method, outdir, readsfile, sitesfile)
