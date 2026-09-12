# Repository guidance

Read [AGENTS.md](AGENTS.md) first. It is the shared entry point for coding and
analysis agents, including Claude Code. The active pattern-model plan and its
source-pinned evidence map are linked there.

Do not infer an active branch, available API or completed test from historical
prose. The previous guidance is preserved at
[adocs/archive/CLAUDE-pre-agent-2026-09-12.md](adocs/archive/CLAUDE-pre-agent-2026-09-12.md)
for history, not as current instructions.

rMVPA remains an R/S3 analysis package. Existing model, dataset, design, global,
regional and searchlight APIs are authoritative; the new agent-state and crossform
adapter protocols are design work until implemented and tested. The numerical
pattern_model core stays in rMVPA. No universal replacement runner is requested.

ITEM continues to delegate trial-covariance numerics to optional fmrilss. Continuous
hrfdecoder materials stay under archive/ and outside active package builds. New
work must preserve those boundaries and avoid monkey-patching or installing into
the user's R library.
