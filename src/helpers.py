""" Helper functions and data types """

import os

__author__ = "Vince Reuter"
__email__ = "vince.reuter@gmail.com"


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
