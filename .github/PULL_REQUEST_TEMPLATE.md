## Summary

<!-- What does this PR do? Link the issue it resolves if applicable: Closes #N -->

## Type of change

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Breaking change (fix or feature that changes existing behaviour)
- [ ] Documentation update

## Changes

<!-- Bullet list of specific changes -->

## Numerical verification (if applicable)

<!-- For solver changes: confirm physical correctness -->
- [ ] Analytical solution matches within expected tolerance
- [ ] No NaN / Inf produced in tested parameter range
- [ ] CFL / stability conditions documented

## Checklist

- [ ] `poetry run pytest tests/ -v` passes
- [ ] `poetry run ruff check src/ tests/` passes
- [ ] `poetry run black --check src/ tests/` passes
- [ ] Docstrings updated for any changed public API
- [ ] CHANGELOG.md entry added under `[Unreleased]`
