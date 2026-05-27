"""Shared pytest configuration for local source-tree imports.

NOTE: This test suite requires ``--import-mode=importlib`` (set in pyproject.toml).
Test subdirectories intentionally omit ``__init__.py`` to avoid shadowing stdlib
modules (e.g. ``tests/io/`` would shadow the built-in ``io`` package with the
default prepend import mode).
"""

import sys
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
