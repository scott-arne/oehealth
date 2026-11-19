# Polymorphic assessment objects with a shared context

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 2)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Polymorphic assessment objects with a shared context.** An abstract
`Assessment` class whose `Assess` takes a `StructureContext`, one common
`AssessmentResult` type, and suites that are themselves assessments.
Rejected: an enum plus central dispatcher (every new check would edit the
enum, dispatcher, serializer, and bindings, which scales poorly for a
growing catalog); a pipeline with finding sinks and injected reference
providers (lifetimes, dependency ordering, and provider interfaces that
five simple checks do not need); a separate result type per assessment
(each needs its own SWIG wrapping and JSON serializer, and suites would
need type erasure anyway). The context is included now because changing a
virtual signature later means editing every assessment.

## What implementation changed

No change.
