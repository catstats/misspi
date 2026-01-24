## Test environments
- Local macOS (R CMD check --as-cran)

## Notes to CRAN
This update addresses the "Additional issues: donttest" reported in CRAN checks.
The default number of parallel workers is now capped when `ncore` is not provided
and respects `_R_CHECK_LIMIT_CORES_`. Examples were updated to run safely during
CRAN checks using single-core execution.
