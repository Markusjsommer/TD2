"""Keep TD2.__version__ aligned with the version declared in setup.py."""
import re
from pathlib import Path

import TD2


def _setup_py_version() -> str:
    setup_text = (Path(__file__).resolve().parent.parent / "setup.py").read_text()
    match = re.search(r'version="([^"]+)"', setup_text)
    assert match, "could not find version= in setup.py"
    return match.group(1)


def test_version_exposed_on_package():
    assert TD2.__version__ == _setup_py_version()


def test_version_is_111_or_newer():
    # Sanity check: the new release line.
    major, minor, *_ = (int(x) for x in TD2.__version__.split("."))
    assert (major, minor) >= (1, 1)
