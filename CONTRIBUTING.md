# Contributing to scop

This guide records the development conventions used in this repository. Keep changes focused, reuse existing helpers, and make behavior visible in the code that owns it.

## Development principles

- Reuse existing helpers and validate at necessary boundaries. Avoid repeating defensive checks for the same invariant.
- Do not add silent fallbacks or hidden behavior switches. Do not add `Sys.getenv()` switches or hidden `options()` switches. If a temporary option change is necessary, save its prior value and restore it with `on.exit()` in the same local scope, and explain why it is needed.
- Do not add unconditional global `options()` or `Sys.setenv()` changes during package load.
- Keep general refactors and feature development in separate pull requests.

## Logging and package startup

- In ordinary business code, use the existing `log_message()` helper for user-facing progress and status messages. Report errors with its error type so they terminate, and report warnings with its warning type so they retain warning semantics.
- `R/scop-package.R` is an explicit startup exception: code there is not required to use `log_message()`. Package startup notices may use `packageStartupMessage()` with `cli` formatting. Do not add a logging proxy, a dependency-availability probe, or multiple fallbacks to make startup messages go through `log_message()`.
- `print.scop_logo()` uses `cat()` for normal print output. Keep that behavior; it is not a business log message.
- `.onAttach()` currently clears `options(scop_env_cache = NULL)`. This is an existing reset of SCOP's process-level Python environment cache. The cache records the environment specification and can include its resolved Python path; `PrepareEnv()` validates and reuses matching cached environment information, while Python-backed helpers may read the cached Python path or module list. Keep this behavior as-is in this documentation change. Do not move it to `.onLoad()`: that would still modify a global option while changing when it is modified. See [R/scop-package.R](R/scop-package.R) and [R/PrepareEnv.R](R/PrepareEnv.R).

## Runtime-optional backends

- Runtime-optional backends such as CellChat and Presto are resolved when their analysis is requested and are not declared package dependencies. The current `Remotes` field in [DESCRIPTION](DESCRIPTION) contains only `thisplot` and `thisutils`, which are imported helpers, not an inventory of optional backends. Do not add a runtime-optional backend to `Imports`, `Depends`, or `Suggests` just to silence an `R CMD check` warning.
- Before executing a runtime-optional backend, use `check_r()` or its existing `<backend>_check_r()` helper. `check_r()` defaults to `install = TRUE`; `verbose = FALSE` only suppresses messages and does not disable installation. Resolve optional backend functions with `get_namespace_fun("<pkg>", "<fun>")` before calling them.
- Runnable examples, tests, and read-only availability checks must not rely on automatic installation. Use `check_r("<repo-or-package>", install = FALSE, verbose = FALSE)` (or a helper with equivalent non-installing behavior), inspect its returned status, and skip or stop clearly when the backend is unavailable. Install needed backends separately as an explicit setup step before running examples or checks.
- For these runtime-optional backends, executable package code and runnable examples must not use `pkg::fun()`, `pkg:::fun()`, `library()`, `require()`, or `requireNamespace()`. This convention does not apply to declared dependencies such as `thisplot` and `thisutils`. Explanatory prose may mention `pkg::fun()`; keep that distinct from executable code.
- When adding an optional backend, centralize its repository/package names, availability check, and function lookup in `<backend>_check_r()` and `<backend>_get_fun()` helpers.
- If roxygen examples change, regenerate the matching `man/*.Rd` files.

## Pull requests and documentation

Keep each pull request scoped to one coherent purpose. Edit source documentation rather than generated website files under `docs/`. Review the final diff for scope and consistency before opening a pull request.
