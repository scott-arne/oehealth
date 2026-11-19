# JSON serialization ships in this round

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 8)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**JSON serialization ships in this round.** nlohmann/json is a private
dependency used only in `src/json.cpp`. Rejected: deferring serialization
(the schema would be retrofitted after reports exist); adding a text
report now (its format would be designed before seeing real results).

## What implementation changed

No change.
