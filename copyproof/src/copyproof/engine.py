"""Find, build and run the k-mer engine (the ``cpk`` binary in engine/)."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

_REPO_ENGINE = Path(__file__).resolve().parents[2] / "engine"


def engine_path(build: bool = True) -> str:
    """Path of the cpk binary: $COPYPROOF_ENGINE, then the source tree, then PATH."""
    env = os.environ.get("COPYPROOF_ENGINE")
    if env:
        if not Path(env).is_file():
            raise FileNotFoundError(f"COPYPROOF_ENGINE={env} does not exist")
        return env
    local = _REPO_ENGINE / "target" / "release" / "cpk"
    if (_REPO_ENGINE / "Cargo.toml").is_file():
        if build and shutil.which("cargo"):
            newest = max(p.stat().st_mtime for p in (_REPO_ENGINE / "src").glob("*.rs"))
            if not local.is_file() or local.stat().st_mtime < newest:
                subprocess.run(["cargo", "build", "--release", "--quiet"], cwd=_REPO_ENGINE, check=True)
        if local.is_file():
            return str(local)
    found = shutil.which("cpk")
    if found:
        return found
    raise FileNotFoundError(
        "the cpk engine was not found: build it with `cargo build --release` in copyproof/engine, "
        "or set COPYPROOF_ENGINE to the binary"
    )


def run_engine(args: list[str], stdin=None, log=None, check: bool = True) -> subprocess.CompletedProcess:
    """Run ``cpk`` with the given arguments. The engine logs progress on stderr."""
    cmd = [engine_path()] + [str(a) for a in args]
    err = log if log is not None else sys.stderr
    return subprocess.run(cmd, stdin=stdin, stderr=err, check=check)
