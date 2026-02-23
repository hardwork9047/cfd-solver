# VIBES.md: Security & Development Integrity

This document outlines the core principles and boundaries for the CFD Solver Tutorial project, ensuring it remains secure, stable, and true to its educational mission.

## 🛡️ Security Practices
- **No Secrets:** Never hardcode API keys, credentials, or personal information. Use environment variables for any external integrations.
- **Dependency Safety:** Only use verified packages from PyPI. Always check `uv.lock` for reproducible and safe environments.
- **Data Integrity:** Protect numerical results. Ensure simulation data is stored in designated directories and not committed unless necessary.

## 💎 Code & Development Integrity
- **Verification:** All code must pass `pytest` and maintain numerical stability (no NaNs or Infs). Refer to `CONTEXT.md` for testing infrastructure.
- **Style:** Adhere strictly to PEP 8 and the established CFD nomenclature (`u`, `v`, `rho`, `mu`). For detailed implementation standards, consult `GEMINI.md`.
- **Transparency:** Document all numerical shortcuts and assumptions.

## 🌟 Professional Demeanor
- **Engineering Excellence:** Write code that is functional, beautiful, and maintainable.
- **Stability First:** Prioritize fixing physics and stability over forcing results.
- **Quality Focus:** Every iteration aims for simulation perfection.
