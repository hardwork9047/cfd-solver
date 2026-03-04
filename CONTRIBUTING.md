# Contribution Guide

## Setup

```bash
# Clone the repository
git clone https://github.com/hardwork9047/cfd-solver.git
cd cfd-solver

# Setup the project with Poetry
poetry install
```

## Development Workflow

### Running Tests

```bash
# Run all tests
poetry run pytest tests/

# Run tests with coverage (all packages)
poetry run pytest tests/ --cov=cfd --cov=cfd_lbm --cov=dem --cov=cfd_dem_lbm --cov-report=html

# Run specific tests
poetry run pytest tests/test_cfd.py::TestPlanePoiseuille -v
```

### Code Quality Checks

```bash
# Check code formatting (black)
poetry run black --check src/ tests/

# Linting (ruff)
poetry run ruff check src/ tests/

# Automatically fix formatting and linting issues
poetry run black src/ tests/
poetry run ruff check --fix src/ tests/
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
   poetry run pytest tests/ -v
   ```

4. **Update Documentation**
   - Update `README.md`
   - Add code comments/docstrings
   - Update `CHANGELOG.md`

5. **Submit a Pull Request**

## Project Structure

```
src/
├── cfd/                     # Navier-Stokes FDM solvers
│   ├── __init__.py
│   ├── poiseuille.py        # Poiseuille flow (Newtonian & Power-Law)
│   ├── cavity.py            # Lid-driven cavity (projection method)
│   ├── cylinder.py          # Flow past a cylinder (Kármán vortex)
│   └── non_newtonian.py     # Non-Newtonian fluid solvers
├── cfd_lbm/                 # Lattice-Boltzmann Method (D2Q9)
│   ├── __init__.py
│   ├── lbm.py               # LBM lid-driven cavity solver
│   └── poiseuille.py        # LBM Poiseuille channel flow
├── dem/                     # Discrete Element Method
│   ├── __init__.py
│   └── particles.py         # DEM particle system (Hertz contact)
└── cfd_dem_lbm/             # Coupled LBM-DEM solver
    ├── __init__.py
    └── lbm_dem.py

tests/
├── test_cfd.py              # Tests for cfd package
├── test_lbm.py              # Tests for cfd_lbm package
└── test_dem.py              # Tests for dem package

examples/                    # Runnable example scripts
.github/
├── workflows/ci.yml         # GitHub Actions CI (tests + lint)
└── ISSUE_TEMPLATE/          # Bug and feature request templates
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

The project uses **Poetry** (`poetry-core`) as its build backend (defined in `pyproject.toml`):
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

### Alternative: Manual Build with Python

If you prefer not to use Poetry for building:

```bash
# Build the package using Python's build module
pip install build
python -m build

# Install the wheel locally for testing
pip install dist/cfd_solver-0.2.0-py3-none-any.whl
```

**Publishing to PyPI with twine:**

```bash
# Install twine
pip install twine

# Upload to PyPI
python -m twine upload dist/*

# You'll be prompted for:
# - Username: __token__
# - Password: <your-pypi-token>
```

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

| Class | Module | Description |
|-------|--------|-------------|
| `PlanePoiseuille` | `cfd` | Poiseuille flow between parallel plates |
| `CircularPoiseuille` | `cfd` | Hagen-Poiseuille flow in a pipe |
| `PowerLawPlanePoiseuille` | `cfd` | Power-law non-Newtonian fluid flow |
| `CavityFlow` | `cfd` | Lid-driven cavity (projection method) |
| `CylinderFlow` | `cfd` | Flow past a cylinder with Kármán vortex shedding |
| `LBM` | `cfd_lbm` | D2Q9 Lattice-Boltzmann cavity solver |
| `PoiseuilleFlow` | `cfd_lbm` | LBM channel (Poiseuille) flow |
| `ParticleSystem` | `dem` | DEM granular particle dynamics |

Refer to [README.md](README.md) for further details.

## Troubleshooting

### If tests fail:

```bash
# Show detailed error messages
poetry run pytest tests/ -vv --tb=long

# Debug a specific test
poetry run pytest tests/test_cfd.py::TestPlanePoiseuille::test_initialization -vv
```

### If you encounter import errors:

```bash
# Rebuild the virtual environment
poetry install --sync
```

## Support

Please use GitHub Issues for questions or bug reports.
