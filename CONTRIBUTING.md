# Contribution Guide

## Setup

```bash
# Clone the repository
git clone https://github.com/hardwork9047/cfd-solver.git
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

This project supports both **Poetry** and **uv** for building and publishing.

### Understanding the Build System

The project uses `hatchling` as its build backend (defined in `pyproject.toml`):
- **Source Distribution (sdist)**: `.tar.gz` file containing source code
- **Wheel**: `.whl` file - a built, ready-to-install package

Both Poetry and uv create these artifacts in the `dist/` directory.

### Method 1: Using Poetry (Recommended for Publishing)

Poetry provides an integrated workflow for building and publishing:

```bash
# Build the package (creates both sdist and wheel)
poetry build

# Output:
# Building cfd-solver (0.2.0)
# Building sdist
# Building wheel
#
# Creates:
# - dist/cfd_solver-0.2.0.tar.gz
# - dist/cfd_solver-0.2.0-py3-none-any.whl
```

**Publishing to PyPI with Poetry:**

```bash
# Option 1: Configure PyPI token once
poetry config pypi-token.pypi <your-pypi-token>
poetry publish

# Option 2: Provide credentials directly
poetry publish --username __token__ --password <your-pypi-token>

# Option 3: Use environment variable
export POETRY_PYPI_TOKEN_PYPI=<your-pypi-token>
poetry publish
```

**Getting a PyPI Token:**
1. Go to https://pypi.org/manage/account/token/
2. Create a new API token
3. Set scope to "Entire account" or specific project
4. Copy the token (starts with `pypi-`)

### Method 2: Using uv

```bash
# Build the package
uv build

# Install the wheel locally for testing
pip install dist/cfd_solver-0.2.0-py3-none-any.whl
```

**Publishing to PyPI with uv and twine:**

```bash
# Install twine
pip install twine

# Upload to PyPI
python -m twine upload dist/*

# You'll be prompted for:
# - Username: __token__
# - Password: <your-pypi-token>
```

### Build Workflow Comparison

| Step | Poetry | uv + twine |
|------|--------|------------|
| Build | `poetry build` | `uv build` |
| Publish | `poetry publish` | `twine upload dist/*` |
| Credentials | Built-in config | Separate tool |
| Speed | Moderate | Fast build |

### Complete Release Workflow

```bash
# 1. Update version in pyproject.toml
# version = "0.3.0"

# 2. Update CHANGELOG.md with new changes

# 3. Commit and tag
git add pyproject.toml CHANGELOG.md
git commit -m "chore: Bump version to 0.3.0"
git tag v0.3.0
git push origin main --tags

# 4. Build and publish
poetry build
poetry publish

# 5. Verify on PyPI
# https://pypi.org/project/cfd-solver/
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
