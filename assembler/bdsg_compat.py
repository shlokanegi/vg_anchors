"""Compatibility checks between the installed libbdsg and the .dist index format.

libbdsg only started recording a version number inside the distance index in
November 2025. A libbdsg built before that reads a newer index with the old
record layout: deserialization succeeds, and the mismatch only surfaces later as
an opaque "traversing the wrong kind of net_handle_t" error from deep inside the
snarl-tree traversal. Catching the stale library up front turns that into
something actionable.
"""

import mmap
import os

# Compiled into libbdsg only once the distance-index version field exists.
_VERSION_CHECK_MARKER = b"up-to-date v"

_REBUILD_INSTRUCTIONS = """Rebuild libbdsg from the pinned submodule:

    git submodule update --init --recursive libbdsg
    pip install -e libbdsg

The pinned commit is the one vg 1.74 builds against, so indexes written by that
vg are readable as-is. If you generate indexes with a different vg, pin libbdsg
to the commit that vg uses (its deps/libbdsg submodule)."""

_STALE_LIBBDSG_MESSAGE = """The installed libbdsg is too old to read this distance index safely.

It predates the distance-index version field, so it cannot report that a .dist
file is in a newer format. It would load the index and then fail while walking
the snarl tree.

{instructions}

Loaded libbdsg: {path}"""

_TRAVERSAL_FAILURE_HINT = """This usually means the .dist index was written by a newer vg than the libbdsg
this tool was built against, so the snarl-tree records are being read with the
wrong layout.

{instructions}

Alternatively, regenerate the index with a vg that matches the installed
libbdsg."""


def _loaded_libbdsg_path():
    """Locate the libbdsg shared object currently mapped into this process.

    Returns
    -------
    str or None
        Path to the mapped libbdsg.so, or None if it cannot be determined.
    """
    try:
        with open("/proc/self/maps", encoding="utf-8") as maps:
            for line in maps:
                fields = line.split()
                if len(fields) < 6:
                    continue
                path = fields[-1]
                if path.startswith("/") and os.path.basename(path).startswith("libbdsg.so"):
                    return path
    except OSError:
        return None
    return None


def _supports_versioned_index(library_path) -> bool:
    """Check whether a libbdsg binary knows about distance-index versions."""
    with open(library_path, "rb") as handle:
        with mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as mapped:
            return mapped.find(_VERSION_CHECK_MARKER) != -1


def check_libbdsg_index_support() -> None:
    """Fail early if the loaded libbdsg cannot detect index version mismatches.

    Skips the check silently when the library cannot be located or read, so this
    never blocks a working setup on an unexpected platform.

    Returns
    -------
    None

    Raises
    ------
    RuntimeError
        If the loaded libbdsg predates the distance-index version field.
    """
    library_path = _loaded_libbdsg_path()
    if library_path is None:
        return

    try:
        supported = _supports_versioned_index(library_path)
    except (OSError, ValueError):
        return

    if not supported:
        raise RuntimeError(
            _STALE_LIBBDSG_MESSAGE.format(
                instructions=_REBUILD_INSTRUCTIONS, path=library_path
            )
        )


def explain_traversal_failure(error):
    """Add the likely cause to a snarl-tree traversal error, when recognised.

    Parameters
    ----------
    error : Exception
        Error raised by SnarlDistanceIndex.traverse_decomposition.

    Returns
    -------
    str or None
        An augmented message, or None if the error is not a known symptom of a
        libbdsg/index mismatch.
    """
    message = str(error)
    if "wrong kind of net_handle_t" not in message:
        return None

    hint = _TRAVERSAL_FAILURE_HINT.format(instructions=_REBUILD_INSTRUCTIONS)
    return f"{message}\n\n{hint}"
