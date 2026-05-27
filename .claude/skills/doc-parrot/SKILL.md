---
name: doc-parrot
description: Validate docstring–implementation alignment for every changed callable in the current branch diff
allowed-tools: Read, Bash, Glob, Grep
---

# Doc-Parrot

You are the doc-parrot agent. Your job is to find every callable that changed in this branch,
extract its docstring prose, and answer one question: **could a developer read this docstring
and use this callable correctly?**

There are two ways a docstring fails that test:
1. It says something **wrong** — the prose contradicts the code.
2. It says something **incomplete** — it's accurate as far as it goes, but omits something a
   caller would reasonably assume, leading them to write code that fails or misbehaves.

The parrot checks both. Do not look at the code until step 5.

## Step 1 + 2: Collect changed callables and extract docstrings

Run this single script to find every callable touched in the branch diff and print its docstring prose.
It keeps lineno as an int throughout to avoid type comparison bugs.

Write the diff to a temp file first — piping directly into `python3 -` while also using a heredoc
conflicts on stdin (the heredoc wins and the diff is lost).

```bash
BASE=$(git merge-base HEAD main)
DIFF_FILE=$(mktemp /tmp/parrot_diff_XXXXX.patch)
git diff "$BASE" HEAD --unified=0 -- '*.py' > "$DIFF_FILE"
python3 << 'PYEOF'
import re, ast

with open("$DIFF_FILE") as f:
    diff_text = f.read()

# Parse diff into per-file sections; collect (filepath, name, lineno) tuples
file_sections = re.split(r'^diff --git ', diff_text, flags=re.MULTILINE)[1:]
candidates = []  # list of (filepath: str, name: str, lineno: int)

for section in file_sections:
    lines = section.splitlines()
    header_match = re.match(r'a/(\S+) b/(\S+)', lines[0])
    if not header_match:
        continue
    filepath = header_match.group(2)
    if not filepath.endswith('.py'):
        continue

    current_new_line = 0
    for line in lines:
        hunk_match = re.match(r'^@@ -\d+(?:,\d+)? \+(\d+)(?:,\d+)? @@', line)
        if hunk_match:
            current_new_line = int(hunk_match.group(1)) - 1
            continue
        if line.startswith('+++') or line.startswith('---'):
            continue
        if not line.startswith('-'):
            current_new_line += 1
        if line.startswith('+') and re.match(r'^\+\s*(def |class )\w', line):
            code_line = line[1:]
            name_match = re.match(r'\s*(def|class)\s+(\w+)', code_line)
            if name_match:
                candidates.append((filepath, name_match.group(2), current_new_line))

# Now extract docstrings via AST — lineno stays int throughout
for filepath, name, lineno in candidates:
    try:
        with open(filepath) as f:
            source = f.read()
        tree = ast.parse(source)
    except Exception:
        continue
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if node.name == name and node.lineno == lineno:
                docstring = ast.get_docstring(node)
                if docstring:
                    print(f"=== {filepath}:{lineno} {name} ===")
                    print(docstring.strip())
                    print()
PYEOF
rm "$DIFF_FILE"
```

Each block beginning with `=== path:line name ===` is one callable's docstring prose.
Callables with no docstring are silently skipped — they are out of scope for the parrot.

## Step 3: Generate descriptions from prose alone

For **each docstring** you collected, do the following **without looking at the implementation**:

> Given only the docstring text below (no method name, no type signatures), write:
> 1. One sentence describing what this callable does.
> 2. A concrete usage example — 3–6 lines of Python showing how a caller would use it.
>
> Do not infer anything from the callable's name. Derive everything from the prose.

Record your output for each callable before moving on.

## Step 4: Adversarial gap-finding

Still working from prose alone — do not read the implementation yet.

For each docstring, ask yourself: **"What would a typical caller be surprised to discover about
this callable that the prose doesn't warn them about?"**

The filter is *surprise*, not *completeness*. The function name and parameter names already
communicate a lot — only flag what those signals don't cover. A function called `assign_robot`
that mutates a dict called `assignments` does not need to document that mutation; the name
makes it obvious. A function called `calculate_distance` that mutates an input would be
surprising and worth flagging.

Look specifically for these classes of gap:

| Gap class | Flag when… |
|---|---|
| **Preconditions** | An input constraint would surprise a caller — not when it's implied by the name or type (e.g. `speed_ms > 0` is not obvious; `non-empty list for a function called filter_X` is) |
| **Edge cases** | Common boundary values (0, empty, `None`, single-element) produce undocumented errors or silent surprises |
| **Error behaviour** | The function can raise but the docstring says nothing about it, or says "raises X" without saying when |
| **Mutation** | The callable modifies an input **and** the function name doesn't already imply that (setters, appenders, and assigners are expected to mutate) |
| **Return ambiguity** | The return description could be read two ways that lead to different calling code |
| **Ordering / idempotency** | Must be called after something else, or calling it twice has different results from calling it once |
| **Side effects** | Writes to disk, logs, mutates shared state, makes a network call — when not implied by the name |

Write down every gap you find — even uncertain ones. You will verify which are real in step 6.
If you find no gaps, write "no gaps identified" and continue.

## Step 5: Run your generated example

For each callable, write the usage example you produced in step 3 to a temp file and run it.
The goal is to check whether the docstring gives a caller enough information to write working code.

```bash
EXAMPLE_FILE=$(mktemp /tmp/parrot_example_XXXXX.py)
cat > "$EXAMPLE_FILE" << 'EOF'
# Insert the import and usage example you generated from prose alone.
# Example:
#   import sys; sys.path.insert(0, '.')
#   from scripts.demo_routing_utils import nearest_slot
#   result = nearest_slot("B04-1", ["B01-1", "B06-1"])
#   print(result)
EOF
python3 "$EXAMPLE_FILE"
rm "$EXAMPLE_FILE"
```

Adjust the import path to match where the callable lives. Record what actually happens:

- **Runs and output matches prose expectations** → the docstring gives a caller enough to work with
- **Runs but output is surprising** → the docstring is technically consistent but misleading about outcomes; flag as **Fix docstring**
- **Errors on import or call** → the docstring described a call signature that doesn't work; flag as **Fix docstring**

## Step 6: Compare against the implementation

Now read the actual implementation of each callable (the code, the tests that cover it).
Bring your description (step 3), your gap list (step 4), and the example run result (step 5).

**Question A — accuracy:** Does the prose say anything wrong? Check the description you derived
against what the code actually does.

**Question B — completeness:** For each gap you flagged in step 4, confirm whether the code
makes it matter. A gap is a real finding only if the code exhibits the behaviour the prose
failed to warn about:
- The undocumented precondition is real and will silently break a caller who violates it
- The undescribed edge case produces surprising output or an undocumented error
- The mutation or side effect actually happens
- The ambiguity resolves to a reading the prose didn't make clear

Assign one judgment per callable:

| Judgment | Meaning |
|---|---|
| **Fix docstring (wrong)** | The prose says something that contradicts the code |
| **Fix docstring (gap)** | The prose is accurate but omits something a caller would need to know |
| **Fix code** | The code diverges from what the docstring promises — the implementation is wrong |
| **OK** | A caller reading only the docstring would use this correctly |

Minor wording differences are OK. Flag things a caller would actually get wrong, not style nits.

## Step 7: Produce the parrot artifact

Output a summary in this format (used by the implementing agent to populate the findings file):

```
## Doc Parrot Results

| Callable | File | Judgment | Notes |
|---|---|---|---|
| foo() | path/to/file.py:42 | OK | — |
| Bar.process() | path/to/other.py:88 | Fix docstring (wrong) | Says "returns a list" but returns a dict |
| Baz.run() | path/to/other.py:120 | Fix docstring (gap) | Mutates input list in-place; not mentioned |

Divergences Found: N
```

If no callables had docstrings, output:
```
## Doc Parrot Results
No docstrings found in changed callables. Nothing to review.
Divergences Found: 0
```

## Notes

- Skip test functions (`test_*`, `Test*`), `__init__`, `__repr__`, `__str__`, `__eq__`, and
  other dunder methods — their contracts are implicit and not worth parroting.
- If a callable is a stub or abstract method with a one-liner docstring like `"""Not implemented."""`,
  skip it.
- The parrot is not a gate. Its output is an artifact the implementing agent acts on with judgment.
  You may decide the divergence is intentional (e.g., the docstring is intentionally simplified
  for external consumers) and record "OK" with a note.
