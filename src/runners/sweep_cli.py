#!/usr/bin/env python3
"""Multi-machine sweep orchestration CLI (invoked via ``bin/sweep``).

Orchestrates manually staged case sweeps across multiple machines using git
as the coordination channel: each machine runs its assigned cases and pushes
lightweight result artifacts to ``sweep_results/<sweep>/<case>/`` so every
machine sees cross-machine progress with a plain ``git pull``.

Stdlib-only on purpose — ``list`` and ``status`` must work on a fresh
checkout before ``poetry install``.

Subcommands:
    list                          show registered sweeps and their cases
    status [<sweep>] [--no-pull]  progress table (done/in-progress/failed/pending)
    run <sweep> <case...|all>     run cases sequentially, harvest + push after each
    harvest <sweep> <case...|all> re-harvest from local results without running
"""

from __future__ import annotations

import argparse
import json
import shutil
import socket
import subprocess
import sys
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFESTS_DIR = REPO_ROOT / "configs" / "lbm_dem" / "sweeps" / "manifests"
RESULTS_ROOT = REPO_ROOT / "src" / "particulate_flow" / "results" / "run_lbm_dem"
SHARED_ROOT = REPO_ROOT / "sweep_results"

# Copied from a completed run directory into sweep_results/<sweep>/<case>/.
# Relative paths; missing entries are skipped silently. Heavy artifacts
# (paraview/, *.png, *.mp4) are intentionally excluded.
HARVEST_FILES = [
    "run_status.json",
    "summary.json",
    "metadata.json",
    "config.json",
    "environment.json",
    "ANALYSIS_REPORT.md",
    "analysis/time_series.csv",
    "analysis/summary.json",
]

STARTED_MARKER = "STARTED.json"
FAILED_MARKER = "FAILED.json"


# ---------------------------------------------------------------------------
# Manifest handling
# ---------------------------------------------------------------------------


def load_manifest(name: str, manifests_dir: Path = MANIFESTS_DIR) -> dict:
    """Load and validate a sweep manifest by name (filename stem)."""
    path = manifests_dir / f"{name}.json"
    if not path.is_file():
        available = ", ".join(m.stem for m in sorted(manifests_dir.glob("*.json")))
        raise SystemExit(
            f"sweep: manifest '{name}' not found in {manifests_dir}"
            + (f" (available: {available})" if available else "")
        )
    manifest = json.loads(path.read_text())
    cases = manifest.get("cases")
    if not isinstance(cases, list) or not cases:
        raise SystemExit(f"sweep: manifest '{name}' has no 'cases' list")
    seen_aliases: set[str] = set()
    for case in cases:
        if "config" not in case:
            raise SystemExit(f"sweep: manifest '{name}' has a case without 'config'")
        alias = case.get("alias")
        if alias:
            if alias in seen_aliases:
                raise SystemExit(f"sweep: manifest '{name}' has duplicate alias '{alias}'")
            seen_aliases.add(alias)
    manifest.setdefault("name", name)
    return manifest


def case_name_for(config_path: Path) -> str:
    """Resolve the result-directory name for a case config.

    Reads the leaf JSON only (no ``extends`` resolution): the result tag must
    be set in the leaf config for the sweep CLI to find its results.
    """
    data = json.loads(config_path.read_text())
    tag = data.get("outputs", {}).get("result_tag") or data.get("name")
    return tag or config_path.stem


def resolve_cases(manifest: dict, selectors: list[str]) -> list[dict]:
    """Map user-supplied selectors (alias, case name, or 'all') to case entries.

    Returns dicts with keys: alias, config (Path), case_name.
    """
    entries = []
    for case in manifest["cases"]:
        config_path = REPO_ROOT / case["config"]
        if not config_path.is_file():
            raise SystemExit(f"sweep: config not found: {case['config']}")
        entries.append(
            {
                "alias": case.get("alias", ""),
                "config": config_path,
                "case_name": case_name_for(config_path),
            }
        )

    if len(selectors) == 1 and selectors[0] == "all":
        return entries

    by_key = {}
    for entry in entries:
        if entry["alias"]:
            by_key[entry["alias"]] = entry
        by_key[entry["case_name"]] = entry

    selected = []
    for sel in selectors:
        if sel not in by_key:
            raise SystemExit(
                f"sweep: unknown case '{sel}' (aliases: "
                + ", ".join(e["alias"] or e["case_name"] for e in entries)
                + ")"
            )
        if by_key[sel] not in selected:
            selected.append(by_key[sel])
    return selected


# ---------------------------------------------------------------------------
# Status
# ---------------------------------------------------------------------------


def case_status(sweep_name: str, case_name: str, shared_root: Path = SHARED_ROOT) -> dict:
    """Classify a case from its shared directory.

    Returns {"status": done|failed|in-progress|pending, "machine": str, "when": str}.
    """
    case_dir = shared_root / sweep_name / case_name
    run_status = case_dir / "run_status.json"
    if run_status.is_file():
        try:
            completed = json.loads(run_status.read_text()).get("status") == "completed"
        except json.JSONDecodeError:
            completed = False
        if completed:
            machine, when = "-", "-"
            harvest = case_dir / "harvest.json"
            if harvest.is_file():
                info = json.loads(harvest.read_text())
                machine = info.get("machine", "-")
                when = info.get("harvested_at", "-")
            return {"status": "done", "machine": machine, "when": when}
    failed = case_dir / FAILED_MARKER
    if failed.is_file():
        info = json.loads(failed.read_text())
        return {
            "status": "failed",
            "machine": info.get("machine", "-"),
            "when": info.get("failed_at", "-"),
        }
    started = case_dir / STARTED_MARKER
    if started.is_file():
        info = json.loads(started.read_text())
        return {
            "status": "in-progress",
            "machine": info.get("machine", "-"),
            "when": info.get("started_at", "-"),
        }
    return {"status": "pending", "machine": "-", "when": "-"}


def _short_when(iso: str) -> str:
    """Shorten an ISO timestamp to 'MM-DD HH:MM' for the status table."""
    try:
        return datetime.fromisoformat(iso).strftime("%m-%d %H:%M")
    except ValueError:
        return iso


def print_status_table(manifest: dict) -> None:
    sweep_name = manifest["name"]
    entries = resolve_cases(manifest, ["all"])
    rows = []
    counts = {"done": 0, "in-progress": 0, "failed": 0, "pending": 0}
    for entry in entries:
        st = case_status(sweep_name, entry["case_name"])
        counts[st["status"]] += 1
        rows.append(
            (
                entry["alias"] or "-",
                entry["case_name"],
                st["status"],
                st["machine"],
                _short_when(st["when"]) if st["when"] != "-" else "-",
            )
        )
    print(f"Sweep: {sweep_name} ({len(rows)} cases)")
    headers = ("ALIAS", "CASE", "STATUS", "MACHINE", "WHEN")
    widths = [max(len(headers[i]), *(len(r[i]) for r in rows)) for i in range(len(headers))]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*headers))
    for row in rows:
        print(fmt.format(*row))
    print(
        f"done {counts['done']} / in-progress {counts['in-progress']}"
        f" / failed {counts['failed']} / pending {counts['pending']}"
    )


# ---------------------------------------------------------------------------
# Git helpers
# ---------------------------------------------------------------------------


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", *args], cwd=REPO_ROOT, check=check, capture_output=True, text=True
    )


def git_pull() -> None:
    """Pull latest; tolerate failure (offline machine) with a warning."""
    result = _git("pull", "--rebase", "--autostash", check=False)
    if result.returncode != 0:
        print(
            "sweep: warning: git pull failed, continuing with local state:\n"
            + result.stderr.strip()
        )


def commit_and_push(case_dir: Path, message: str, push: bool = True) -> None:
    """Stage one case directory, commit if anything changed, then push with retry."""
    _git("add", "--", str(case_dir.relative_to(REPO_ROOT)))
    staged = _git("diff", "--cached", "--quiet", check=False)
    if staged.returncode == 0:
        return  # nothing to commit
    _git("commit", "-m", message)
    if not push:
        print("sweep: commit created locally (push skipped)")
        return
    for attempt in range(3):
        result = _git("push", check=False)
        if result.returncode == 0:
            return
        pull = _git("pull", "--rebase", "--autostash", check=False)
        if pull.returncode != 0:
            _git("rebase", "--abort", check=False)
            break
    print(
        "sweep: warning: push failed; the commit stays local and the next"
        " run/harvest will retry the push"
    )


# ---------------------------------------------------------------------------
# Harvest
# ---------------------------------------------------------------------------


def find_completed_run_dir(
    case_name: str, since: float | None = None, results_root: Path = RESULTS_ROOT
) -> Path | None:
    """Newest local run dir for a case whose run_status.json says completed.

    Searches ``dim3/<case>/run_*`` (3D) and ``<case>/run_*`` (2D). ``since``
    is a POSIX timestamp; run dirs whose run_status.json is older are ignored
    (used to pick the run started by this invocation, not a stale one).
    """
    candidates = []
    for parent in (results_root / "dim3" / case_name, results_root / case_name):
        if not parent.is_dir():
            continue
        for run_dir in parent.glob("run_*"):
            status_file = run_dir / "run_status.json"
            if not status_file.is_file():
                continue
            try:
                if json.loads(status_file.read_text()).get("status") != "completed":
                    continue
            except json.JSONDecodeError:
                continue
            if since is not None and status_file.stat().st_mtime < since:
                continue
            candidates.append(run_dir)
    return max(candidates, key=lambda p: p.name) if candidates else None


def harvest_case(
    sweep_name: str, case_name: str, run_dir: Path, shared_root: Path = SHARED_ROOT
) -> Path:
    """Copy lightweight artifacts from ``run_dir`` into the shared case dir."""
    case_dir = shared_root / sweep_name / case_name
    case_dir.mkdir(parents=True, exist_ok=True)
    for rel in HARVEST_FILES:
        src = run_dir / rel
        if not src.is_file():
            continue
        dst = case_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    (case_dir / "harvest.json").write_text(
        json.dumps(
            {
                "machine": socket.gethostname(),
                "source_run_dir": str(run_dir.relative_to(REPO_ROOT)),
                "harvested_at": datetime.now().isoformat(timespec="seconds"),
            },
            indent=2,
        )
        + "\n"
    )
    for marker in (STARTED_MARKER, FAILED_MARKER):
        (case_dir / marker).unlink(missing_ok=True)
    return case_dir


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------


def cmd_list(_args: argparse.Namespace) -> int:
    manifests = sorted(MANIFESTS_DIR.glob("*.json"))
    if not manifests:
        print(f"sweep: no manifests in {MANIFESTS_DIR.relative_to(REPO_ROOT)}")
        return 1
    for path in manifests:
        manifest = load_manifest(path.stem)
        print(f"{manifest['name']}: {manifest.get('description', '')}")
        for case in manifest["cases"]:
            config_path = REPO_ROOT / case["config"]
            alias = case.get("alias", "-")
            if config_path.is_file():
                name = case_name_for(config_path)
                tag_set = bool(
                    json.loads(config_path.read_text()).get("outputs", {}).get("result_tag")
                )
                warn = "" if tag_set else "  [warn: no outputs.result_tag in leaf config]"
                print(f"  {alias:<6} {name}{warn}")
            else:
                print(f"  {alias:<6} [missing config: {case['config']}]")
        print()
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    if not args.no_pull:
        git_pull()
    names = [args.sweep] if args.sweep else sorted(m.stem for m in MANIFESTS_DIR.glob("*.json"))
    if not names:
        print("sweep: no manifests registered")
        return 1
    for name in names:
        print_status_table(load_manifest(name))
        print()
    return 0


def _ensure_poetry_env() -> None:
    print("sweep: poetry install (use --no-install to skip)")
    subprocess.run(["poetry", "install"], cwd=REPO_ROOT, check=True)


def cmd_run(args: argparse.Namespace) -> int:
    manifest = load_manifest(args.sweep)
    cases = resolve_cases(manifest, args.cases)
    sweep_name = manifest["name"]

    if args.dry_run:
        for entry in cases:
            st = case_status(sweep_name, entry["case_name"])
            action = "SKIP (completed)" if st["status"] == "done" and not args.force else "RUN"
            config_rel = entry["config"].relative_to(REPO_ROOT)
            cmd = f"poetry run python src/runners/run_lbm_dem.py --config {config_rel}"
            print(f"{action:<18} {entry['case_name']}\n{'':<18} {cmd}")
        return 0

    git_pull()
    if not args.no_install:
        _ensure_poetry_env()

    failures = 0
    for entry in cases:
        case_name = entry["case_name"]
        st = case_status(sweep_name, case_name)
        if st["status"] == "done" and not args.force:
            print(
                f"SKIP {case_name} (completed on {st['machine']} at {st['when']};"
                " use --force to rerun)"
            )
            continue
        if st["status"] == "in-progress":
            print(
                f"sweep: note: {case_name} is marked in-progress on {st['machine']}"
                f" (started {st['when']}) — running anyway per manual assignment"
            )

        case_dir = SHARED_ROOT / sweep_name / case_name
        case_dir.mkdir(parents=True, exist_ok=True)
        (case_dir / STARTED_MARKER).write_text(
            json.dumps(
                {
                    "machine": socket.gethostname(),
                    "started_at": datetime.now().isoformat(timespec="seconds"),
                    "config": str(entry["config"].relative_to(REPO_ROOT)),
                },
                indent=2,
            )
            + "\n"
        )
        commit_and_push(
            case_dir,
            f"sweep({sweep_name}): start {case_name} on {socket.gethostname()}",
            push=not args.no_push,
        )

        started_at = datetime.now().timestamp()
        print(f"RUN  {case_name}  ({entry['config'].relative_to(REPO_ROOT)})")
        result = subprocess.run(
            [
                "poetry",
                "run",
                "python",
                "src/runners/run_lbm_dem.py",
                "--config",
                str(entry["config"].relative_to(REPO_ROOT)),
            ],
            cwd=REPO_ROOT,
        )

        if result.returncode == 0:
            run_dir = find_completed_run_dir(case_name, since=started_at)
            if run_dir is None:
                print(
                    f"sweep: warning: {case_name} exited 0 but no completed run dir"
                    " found; not harvesting"
                )
                failures += 1
            else:
                harvest_case(sweep_name, case_name, run_dir)
                commit_and_push(
                    case_dir,
                    f"sweep({sweep_name}): complete {case_name} on {socket.gethostname()}",
                    push=not args.no_push,
                )
                print(f"DONE {case_name}  (harvested from {run_dir.relative_to(REPO_ROOT)})")
                continue
        else:
            failures += 1
            (case_dir / FAILED_MARKER).write_text(
                json.dumps(
                    {
                        "machine": socket.gethostname(),
                        "failed_at": datetime.now().isoformat(timespec="seconds"),
                        "exit_code": result.returncode,
                    },
                    indent=2,
                )
                + "\n"
            )
            (case_dir / STARTED_MARKER).unlink(missing_ok=True)
            commit_and_push(
                case_dir,
                f"sweep({sweep_name}): FAILED {case_name} on {socket.gethostname()}",
                push=not args.no_push,
            )
            print(f"FAIL {case_name} (exit {result.returncode})")

        if not args.keep_going:
            print("sweep: stopping after failure (use --keep-going to continue)")
            break

    return 1 if failures else 0


def cmd_harvest(args: argparse.Namespace) -> int:
    manifest = load_manifest(args.sweep)
    cases = resolve_cases(manifest, args.cases)
    sweep_name = manifest["name"]
    harvested = 0
    for entry in cases:
        case_name = entry["case_name"]
        run_dir = find_completed_run_dir(case_name)
        if run_dir is None:
            print(f"skip {case_name}: no completed local run found")
            continue
        case_dir = harvest_case(sweep_name, case_name, run_dir)
        commit_and_push(
            case_dir,
            f"sweep({sweep_name}): harvest {case_name} from {socket.gethostname()}",
            push=not args.no_push,
        )
        print(f"harvested {case_name} from {run_dir.relative_to(REPO_ROOT)}")
        harvested += 1
    print(f"sweep: harvested {harvested}/{len(cases)} cases")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="sweep", description="Multi-machine sweep orchestration")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("list", help="show registered sweeps and their cases")

    p_status = sub.add_parser("status", help="show cross-machine progress")
    p_status.add_argument("sweep", nargs="?", help="sweep name (default: all)")
    p_status.add_argument("--no-pull", action="store_true", help="skip git pull")

    p_run = sub.add_parser("run", help="run cases sequentially, harvest + push after each")
    p_run.add_argument("sweep", help="sweep manifest name")
    p_run.add_argument("cases", nargs="+", help="case aliases/names, or 'all'")
    p_run.add_argument("--force", action="store_true", help="rerun completed cases")
    p_run.add_argument("--dry-run", action="store_true", help="print actions only")
    p_run.add_argument("--no-push", action="store_true", help="commit locally, do not push")
    p_run.add_argument("--no-install", action="store_true", help="skip poetry install")
    p_run.add_argument("--keep-going", action="store_true", help="continue past case failures")

    p_harvest = sub.add_parser("harvest", help="re-harvest completed local results without running")
    p_harvest.add_argument("sweep", help="sweep manifest name")
    p_harvest.add_argument("cases", nargs="+", help="case aliases/names, or 'all'")
    p_harvest.add_argument("--no-push", action="store_true", help="commit locally, do not push")

    args = parser.parse_args(argv)
    handler = {
        "list": cmd_list,
        "status": cmd_status,
        "run": cmd_run,
        "harvest": cmd_harvest,
    }[args.command]
    return handler(args)


if __name__ == "__main__":
    sys.exit(main())
