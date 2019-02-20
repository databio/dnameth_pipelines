""" Helper functions and data types """

import os
import re
import sys
if sys.version_info < (3, 3):
    from collections import Sequence
else:
    from collections.abc import Sequence


__author__ = "Vince Reuter"
__email__ = "vince.reuter@gmail.com"


class EpilogTarget(Sequence):
    """
    Representation of the main epilog targets.

    This class offer flexibility with respect to presence/absence of
    epialleles file and provides named access to the targets while mixing in
    Sequence behavior for use more like a list.
    """

    def __init__(self, single_sites_file, epis_file=None):
        self._ss = single_sites_file
        self._epi = epis_file
        self._files = [epis_file, single_sites_file] if epis_file else [single_sites_file]

    @property
    def single_sites_file(self):
        return self._ss

    @property
    def epis_file(self):
        return self._epi

    def __getitem__(self, item):
        return self._files[item]

    def __len__(self):
        return len(self._files)


class FolderContext(object):
    """ Context manager for temporarily changing directory. """

    def __init__(self, folder):
        """
        Store the previous working path to restore upon exit.

        :param str folder: Path to set as new working directory
        """
        if not os.path.isdir(folder):
            raise ValueError(
                "Requested temp entry to non-folder: {}".format(folder))
        self._prevdir = os.getcwd()
        self._currdir = folder

    def __enter__(self):
        """ Make the working directory switch. """
        os.chdir(self._currdir)

    def __exit__(self, exc_type, exc_val, exc_tb):
        """ Switch back to the previous working directory. """
        if not os.path.isdir(self._prevdir):
            raise RuntimeError("Return path is no longer a directory: {}".
                               format(self._prevdir))
        os.chdir(self._prevdir)


class ProgSpec(object):
    """ Resource specification for a downstream analysis program. """

    def __init__(self, jar, memory, cores=1):
        """
        A JAR path, memory text, and number of processors specifies program resources.

        :param str jar: Path to JAR in which to find program
        :param str memory: Memory specification; include units and be aware
            of your cluster's / computer's units
        :param int cores: Number of cores to use
        """

        self.memory = memory

        if not os.path.isfile(os.path.expanduser(os.path.expandvars(jar))):
            raise MissingEpilogError("Path to JAR isn't a file: {}".format(jar))
        self.jar = jar

        def check_cores(n):
            try:
                n = int(n)
            except:
                return False
            return (n > 0) and n

        cores = check_cores(cores)
        if not cores:
            raise ValueError("Invalid cores specification ({}); provide a nonnegative integer".format(cores))
        self.cores = cores

    def get_command_base(self, pkg_prog_pair=None):
        """
        Get the base for a command assocaited with this program specification.

        :param str, str pkg_prog_pair: Pair of package name and program name; no specification implies main JAR program
        :return str: Base of command associated with this program specification
        """
        base = "java -Xmx{}".format(self.memory)
        if pkg_prog_pair is None:
            return "{b} -jar {j}".format(b=base, j=self.jar)
        else:
            try:
                pkg, prog = pkg_prog_pair
            except (TypeError, ValueError):
                raise ValueError(
                    "If specifying program, provide a 2-tuple of package name and program name; got {} ({})".
                    format(pkg_prog_pair, type(pkg_prog_pair)))
            return "{b} -cp {j} {p}".format(b=base, j=self.jar, p=self._get_program(pkg, prog))

    @staticmethod
    def _get_program(pkg, prog):
        """ Get fully qualified classpath for program to run. """
        return ".".join(["episcall", pkg, prog])


def get_dedup_bismark_cmd(paired, infile,
    outdir=None, prog="deduplicate_bismark", compress=True):
    """
    Create command with which to run deduplication with bismark.

    :param bool paired: Whether the reads are from a paired-end protocol.
    :param str infile: Path to aligned reads file to deduplicate.
    :param str | NoneType outdir: Path to output folder; if unspecified,
        the parent of the infile will be used
    :param str prog: Path to program to run, or name of program on PATH
    :param bool compress: Whether to compress the SAM output (into BAM)
    :return str, str: Command with which to run bismark deduplication, and the
        output file that the command should create once executed
    """
    read_end_type = "paired" if paired else "single"
    if not os.path.isfile(infile):
        raise Exception("Input to {} is not a file: {}".format(prog, infile))
    cmd = "{} --{} {}".format(prog, read_end_type, infile)
    outfile = re.sub(r'.bam$', '.deduplicated.bam', infile)
    outdir = outdir or os.path.dirname(outfile)
    if outdir:
        if not os.path.isdir(outdir):
            raise Exception("Output folder for {} isn't a directory: {}".format(prog, outdir))
        cmd += " --output_dir {}".format(outdir)
    if compress:
        cmd += " --bam"
    return cmd, outfile


def missing_targets(targets, good=lambda f: os.path.isfile(f)):
    """
    Find missing target(s) of a command.

    :param str | list[str] t: Single target path, or collection of them
    :param callable(str) -> bool good: Predicate to evaluate on a path; assumption
        is that each path is a file and therefore default predicate is
        path's existence as a file.
    :return list[str]: Each target that isn't a file
    """
    import sys
    if sys.version_info < (3, 3):
        from collections import Iterable
    else:
        from collections.abc import Iterable
    if isinstance(targets, str):
        targets = [targets]
    elif not isinstance(targets, Iterable):
        raise TypeError(
            "Target(s) to check should be a collection (or maybe string); got {} ({})".format(targets, type(targets)))
    return [p for p in targets if not good(p)]


class MissingEpilogError(Exception):
    """ Exception for when a program specification's JAR is not a file. """
    pass
