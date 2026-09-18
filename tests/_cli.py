"""Run a SpatialAtlas command-line script in-process.

Importing and calling ``main()`` keeps the scripts inside the coverage report,
which a subprocess invocation would not.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

SOFTWARE_DIR = Path(__file__).resolve().parents[1] / "software"


def load_script(stem: str):
    module_name = f"{stem}_under_test"
    if module_name in sys.modules:
        return sys.modules[module_name]
    spec = importlib.util.spec_from_file_location(
        module_name, SOFTWARE_DIR / f"{stem}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def run_main(stem: str, argv: list[str]):
    """Call ``main()`` with the given arguments; return (exit_code, module)."""
    module = load_script(stem)
    saved = sys.argv
    sys.argv = [f"{stem}.py", *argv]
    try:
        module.main()
        return 0, module
    except SystemExit as exc:
        code = exc.code
        if code is None:
            return 0, module
        return (1 if isinstance(code, str) else int(code)), module
    finally:
        sys.argv = saved
