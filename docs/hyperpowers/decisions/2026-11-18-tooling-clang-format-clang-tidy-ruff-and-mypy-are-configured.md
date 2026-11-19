# Tooling: clang-format, clang-tidy, ruff, and mypy are configured from the start

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 14)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Tooling: clang-format, clang-tidy, ruff, and mypy are configured from
the start.** Rejected: relying on tests alone.

## What implementation changed

No change. `.clang-format`, `.clang-tidy` and the `[tool.ruff]` and `[tool.mypy]` sections of `pyproject.toml` were added. `.clang-tidy` disables a few checks that conflict with project conventions, and each is commented in the file.
