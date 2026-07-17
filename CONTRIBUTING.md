# Contributing

Thank you for considering a contribution. This repository is intended to
support transparent, method-level reproduction of manuscript analyses.

## Reporting Issues

Please open an issue for:

- script errors or unclear workflow instructions;
- missing dependency information;
- inconsistent file names, labels or method descriptions;
- questions about reproducing a documented workflow.

Do not post precise cave roost coordinates, sensitive locality details,
private sampling-site information or controlled-access data in public issues.
For sensitive data-access questions, contact Dr. Zhiqiang Wu at
zqwu_lab@163.com.

## Pull Requests

Before submitting a pull request:

1. Keep changes scoped to one workflow or documentation topic when possible.
2. Preserve existing input/output conventions unless the change explicitly
   updates a workflow contract.
3. Do not add sensitive locality data or unreleased raw data.
4. Run the affected script or document why it could not be run.
5. Update the relevant README if behavior, inputs, outputs or method wording
   changes.

## Reproducibility Expectations

Analysis changes should include enough context for a reader to understand:

- which input files are used;
- which script produced each changed output;
- whether numerical results changed;
- whether downstream workflows need to be rerun.

Documentation-only changes should avoid altering generated analysis products.
