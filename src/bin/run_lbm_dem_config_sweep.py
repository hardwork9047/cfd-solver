#!/usr/bin/env python3
"""Run LBM-DEM parameter sweeps from a JSON configuration file."""

from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
RUNNER = REPO_ROOT / "src" / "demos" / "run_lbm_dem.py"
STATUS_HEADER = ["case", "status", "exit_code", "log_file", "command"]
NEGATABLE_OPTIONS = {
    "flow-control",
    "rolling-friction",
    "paraview-output",
}


def _option_name(name: str) -> str:
    return "--" + name.replace("_", "-")


def _append_option(argv: list[str], name: str, value: Any) -> None:
    option = _option_name(name)
    option_key = option.removeprefix("--")
    if isinstance(value, bool):
        if value:
            argv.append(option)
        elif option_key in NEGATABLE_OPTIONS:
            argv.append("--no-" + option_key)
        return
    if value is None:
        return
    if isinstance(value, list):
        if name.replace("_", "-") == "cylinder-spec" and value and isinstance(value[0], list):
            for spec in value:
                argv.append(option)
                argv.extend(str(item) for item in spec)
            return
        for item in value:
            argv.append(option)
            if isinstance(item, list):
                argv.extend(str(part) for part in item)
            else:
                argv.append(str(item))
        return
    argv.extend([option, str(value)])


def _args_from_mapping(mapping: dict[str, Any]) -> list[str]:
    argv: list[str] = []
    for name, value in mapping.items():
        _append_option(argv, name, value)
    return argv


def _load_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        config = json.load(handle)
    if not isinstance(config, dict):
        raise ValueError("configuration root must be a JSON object")
    if "cases" in config and not isinstance(config.get("cases"), list):
        raise ValueError("cases must be a list")
    if "factors" in config and not isinstance(config.get("factors"), dict):
        raise ValueError("factors must be an object")
    if "cases" not in config and "factors" not in config:
        raise ValueError("configuration must contain cases or factors")
    return config


def _factor_cases(config: dict[str, Any]) -> list[dict[str, Any]]:
    factors = config.get("factors", {})
    if not factors:
        return []
    factor_names = list(factors.keys())
    factor_values: list[list[Any]] = []
    for name in factor_names:
        values = factors[name]
        if not isinstance(values, list) or not values:
            raise ValueError(f"factor {name!r} must be a non-empty list")
        factor_values.append(values)

    cases: list[dict[str, Any]] = []
    for combo in itertools.product(*factor_values):
        args = dict(zip(factor_names, combo))
        label_parts = []
        for name, value in args.items():
            safe_value = str(value).replace("/", "_").replace(" ", "_").replace(".", "p")
            label_parts.append(f"{name}_{safe_value}")
        cases.append({"name": "__".join(label_parts), "args": args})
    return cases


def _all_cases(config: dict[str, Any]) -> list[dict[str, Any]]:
    cases = list(config.get("cases", []))
    cases.extend(_factor_cases(config))
    return cases


def _case_command(config: dict[str, Any], case: dict[str, Any]) -> list[str]:
    defaults = config.get("defaults", {})
    case_args = case.get("args", {})
    if not isinstance(defaults, dict) or not isinstance(case_args, dict):
        raise ValueError("defaults and each case args must be JSON objects")
    merged = {**defaults, **case_args}
    name = str(case.get("name", "case"))
    merged.setdefault("result_tag", name)
    return [sys.executable, "-u", str(RUNNER), *_args_from_mapping(merged)]


def _write_status_row(path: Path, row: dict[str, str]) -> None:
    exists = path.exists()
    with path.open("a", encoding="utf-8") as handle:
        if not exists:
            handle.write("\t".join(STATUS_HEADER) + "\n")
        handle.write("\t".join(row[column] for column in STATUS_HEADER) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run LBM-DEM sweeps from JSON config")
    parser.add_argument("config", type=Path, help="JSON sweep configuration")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running them")
    parser.add_argument("--stop-on-failure", action="store_true", help="Stop at first failed case")
    args = parser.parse_args()

    config = _load_config(args.config)
    run_id = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    log_root = REPO_ROOT / "src" / "results" / "run_lbm_dem" / "config_sweeps" / run_id
    log_root.mkdir(parents=True, exist_ok=True)
    status_path = log_root / "status.tsv"

    cases = _all_cases(config)
    if not cases:
        raise ValueError("configuration expanded to zero cases")

    failures = 0
    for idx, case in enumerate(cases, start=1):
        if not isinstance(case, dict):
            raise ValueError("each case must be a JSON object")
        name = str(case.get("name", f"case_{idx:03d}"))
        command = _case_command(config, case)
        log_file = log_root / f"{idx:03d}_{name}.log"
        quoted = " ".join(subprocess.list2cmdline([part]) for part in command)
        print(f"[{idx}/{len(cases)}] {name}")
        print(f"  {quoted}")
        if args.dry_run:
            continue
        with log_file.open("w", encoding="utf-8") as handle:
            handle.write(f"Command: {quoted}\n")
            completed = subprocess.run(command, cwd=REPO_ROOT, stdout=handle, stderr=subprocess.STDOUT)
        status = "ok" if completed.returncode == 0 else "failed"
        if completed.returncode != 0:
            failures += 1
        _write_status_row(
            status_path,
            {
                "case": name,
                "status": status,
                "exit_code": str(completed.returncode),
                "log_file": str(log_file),
                "command": quoted,
            },
        )
        if completed.returncode != 0 and args.stop_on_failure:
            break

    if not args.dry_run:
        print(f"Status: {status_path}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
