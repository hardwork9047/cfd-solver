# Verification Programs

This directory contains manual verification and benchmark-comparison programs.
They are intentionally separate from the fast pytest suite because several
scripts generate figures, write result files, or run longer numerical studies.

Run the consolidated historical comparison with:

```bash
poetry run python tests/verification/compare_historical_results.py
```

For a slower lid-driven-cavity comparison against Ghia et al. (1982), add:

```bash
poetry run python tests/verification/compare_historical_results.py --include-lbm --lbm-steps 5000
```
