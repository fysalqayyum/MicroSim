# Contributing

Contributions are welcome when they are scoped, reproducible, and honest about
what was tested.

## Before changing a solver

1. Identify the exact active host and device code paths.
2. Record a small reference case and its source/executable checksums.
3. Change one concern at a time.
4. Add a regression check that would have caught the original defect.
5. Re-run the relevant qualification sequence.

An input value printed at startup is not proof that the selected kernel uses
it.

## Branch and commit practice

- Rebase a clean feature branch on the current upstream baseline.
- Keep upstream synchronization separate from fork-specific repairs.
- Use focused commit messages that describe behavior, not the editing process.
- Do not commit generated binaries, full simulation results, scheduler logs,
  core dumps, local environments, or editor state.
- Preserve original authorship and GPL-3.0 notices.

## HPC changes

Every HPC result used to justify a change should record:

- compiler, MPI, accelerator runtime, and library versions;
- GPU/CPU model and rank/thread topology;
- build command;
- executable and source hashes;
- exact inputs;
- scheduler state, exit code, and wall time;
- finite-field and equivalence diagnostics.

Performance claims require an equivalent numerical result, not merely a faster
exit.

## Thermodynamic data

Do not submit licensed or access-controlled TDB files. A contribution may
include:

- code that accepts a user-supplied database;
- documented provenance fields;
- generated tables only when redistribution terms are clear;
- small synthetic test data with explicit units.

Do not infer mobility or diffusivity from a Gibbs-energy database.

## Public-tree safety

Before committing:

```bash
git diff --check
scripts/pre-publish-check.sh upstream/main
```

Review the complete outgoing patch:

```bash
git diff --stat upstream/main...HEAD
git diff upstream/main...HEAD
```

The automated scan is a guardrail, not a substitute for human review.
