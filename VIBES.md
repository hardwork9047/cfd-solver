# VIBES.md: Security & Development Integrity

This document outlines the core principles and boundaries for the CFD Solver Tutorial project, ensuring it remains secure, stable, and true to its educational mission.

## 🛡️ Security Practices
- **No Secrets:** Never hardcode API keys, credentials, or personal information. For GitHub operations (like `git push`), use environment variables to store your Personal Access Token (PAT). For example, set the `GITHUB_HARDWORK9047` environment variable with your PAT. You can then configure Git to use this variable for authentication, for instance, by setting a global Git configuration like `git config --global http.https://github.com/.extraheader 'AUTHORIZATION: token $GITHUB_HARDWORK9047'`.
- **Dependency Safety:** Only use verified packages from PyPI. Always check `uv.lock` for reproducible and safe environments.
- **Data Integrity:** Protect numerical results. Ensure simulation data is stored in designated directories and not committed unless necessary.
- **Git commands:** All git commands are allowed via Gemini. (user:hardwork9047, password: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)

## 💎 Code & Development Integrity
- **Verification:** All code must pass `pytest` and maintain numerical stability (no NaNs or Infs). Refer to `CONTEXT.md` for testing infrastructure.
- **Style:** Adhere strictly to PEP 8 and the established CFD nomenclature (`u`, `v`, `rho`, `mu`). For detailed implementation standards, consult `GEMINI.md`.
- **Transparency:** Document all numerical shortcuts and assumptions.

## 🌟 Professional Demeanor
- **Engineering Excellence:** Write code that is functional, beautiful, and maintainable.
- **Stability First:** Prioritize fixing physics and stability over forcing results.
- **Quality Focus:** Every iteration aims for simulation perfection.
