# Security Policy

## Supported Versions

| Version | Supported |
|---------|-----------|
| 0.3.x   | ✅ Yes    |
| 0.2.x   | ✅ Yes    |
| < 0.2   | ❌ No     |

## Reporting a Vulnerability

Please **do not** report security vulnerabilities through public GitHub Issues.

Instead, use one of the following private channels:

1. **GitHub Security Advisories** (preferred):
   Go to the [Security tab](https://github.com/hardwork9047/cfd-solver/security/advisories/new)
   and click "Report a vulnerability".

2. **Email**: Contact the maintainer directly via the email listed on the
   [GitHub profile](https://github.com/hardwork9047).

### What to include

- A clear description of the vulnerability
- Steps to reproduce (minimal reproducible example if possible)
- The potential impact
- Any suggested mitigations

### Response timeline

- **Acknowledgement**: Within 48 hours
- **Initial assessment**: Within 7 days
- **Fix or mitigation**: Within 30 days for confirmed vulnerabilities

We will credit reporters in the release notes unless you prefer to remain anonymous.

## Security Best Practices for Users

- Never hardcode API keys, tokens, or credentials in code that uses this library.
- Pin dependency versions in production environments (`poetry.lock`).
- Validate all external inputs before passing them to solver parameters.
- Simulation results may be sensitive; protect output files appropriately.
