# Contribution Guide

## Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/cfd-solver.git
cd cfd-solver

# Setup the project with uv
uv sync --all-extras
```

## Development Workflow

### Running Tests

```bash
# Run all tests
uv run pytest tests/

# Run tests with coverage
uv run pytest tests/ --cov=cfd --cov-report=html

# Run specific tests
uv run pytest tests/test_cfd.py::TestPlanePoiseuille -v
```

### Code Quality Checks

```bash
# Check code formatting (black)
uv run black --check src/cfd tests/

# Linting (ruff)
uv run ruff check src/cfd tests/

# Automatically fix formatting and linting issues
uv run black src/cfd tests/
uv run ruff check --fix src/cfd tests/
```

## Adding New Features

1. **Create a Feature Branch**
   ```bash
   git checkout -b feature/new-feature
   ```

2. **Implement Your Code**
   - Add new modules to `src/cfd/`
   - Add corresponding tests in `tests/`

3. **Verify with Tests**
   ```bash
   uv run pytest tests/ -v
   ```

4. **Update Documentation**
   - Update `README.md`
   - Add code comments/docstrings
   - Update `CHANGELOG.md`

5. **Submit a Pull Request**

## Project Structure

```
src/cfd/
├── __init__.py          # Package initialization and exports
├── poiseuille.py        # Poiseuille flow implementations
├── cavity.py            # Lid-driven cavity flow implementation
└── non_newtonian.py     # Non-Newtonian fluid solvers

tests/
├── test_cfd.py          # Unit tests
└── conftest.py          # pytest configurations (optional)
```

## Release Procedure

1. Update `CHANGELOG.md`
2. Update the version in `pyproject.toml`
3. Create a git tag: `git tag v0.x.x`
4. Push to GitHub
5. Publish to PyPI (via CI/CD)

## Building and Publishing the Package

### Local Build

```bash
# Create the wheel
uv build

# Install the wheel locally
pip install dist/cfd_solver-0.2.0-py3-none-any.whl
```

### Publishing to PyPI

```bash
# Build the package
uv build

# Upload to PyPI using twine (ensure credentials are set in ~/.pypirc)
python -m twine upload dist/*
```

## Documentation

For detailed API documentation, please refer to the docstrings in each module.

### Core Classes

- `PlanePoiseuille`: Poiseuille flow between parallel plates.
- `CircularPoiseuille`: Hagen-Poiseuille flow in a pipe.
- `PowerLawPlanePoiseuille`: Flow of Power-Law fluids.
- `CavityFlow`: Lid-driven cavity flow simulation.

Refer to [README.md](README.md) for further details.

## Troubleshooting

### If tests fail:

```bash
# Show detailed error messages
uv run pytest tests/ -vv --tb=long

# Debug a specific test
uv run pytest tests/test_cfd.py::TestPlanePoiseuille::test_initialization -vv
```

### If you encounter import errors:

```bash
# Rebuild the virtual environment
uv sync --force-reinstall
```

## Support

Please use GitHub Issues for questions or bug reports.
