#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

echo "Running automated vibe checks..."
echo "-------------------------------"

# Check Ruff (linting and style)
echo "Running Ruff check..."
poetry run ruff check .
if [ $? -ne 0 ]; then
    echo "Ruff check found issues. Please review the output above."
    # In a CI environment, you might want to exit with an error code here:
    # exit 1
fi
echo "Ruff check completed."
echo "-------------------------------"

# Check Black (formatting)
echo "Running Black check..."
poetry run black --check .
if [ $? -ne 0 ]; then
    echo "Black check found formatting issues. Please run 'poetry run black .' to format."
    # In a CI environment, you might want to exit with an error code here:
    # exit 1
fi
echo "Black check completed."
echo "-------------------------------"

# Run Pytest (verification)
echo "Running Pytest..."
poetry run pytest
if [ $? -ne 0 ]; then
    echo "Pytest failed. Please review the test failures."
    # In a CI environment, you might want to exit with an error code here:
    # exit 1
fi
echo "Pytest completed."
echo "-------------------------------"

echo "Automated vibe checks finished. Please review the output above."
echo "To apply fixes, run 'poetry run ruff --fix .' and 'poetry run black .' separately, then commit the changes."

exit 0
