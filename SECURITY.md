# Data and credential safety

Do not open an issue containing participant-level data, restricted SWAN files,
credentials, access tokens, personal contact details, or local configuration.

Before committing, run:

```powershell
Rscript --vanilla tests/audit_release.R
```

The audit rejects common data formats, manuscript formats, generated outputs,
credential patterns, and absolute local paths.

