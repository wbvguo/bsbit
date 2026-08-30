# Fixture maintenance tools

These tools reproduce or independently validate small stable fixtures.
They are not product binaries and are not invoked by a normal Cargo build.
Generated large data or one-run output must be written under ignored
`workspace/`.

The maintained tests use only the Python standard library.

Run every maintained policy/tool test from the repository root:

```sh
python3 -m unittest discover -s tests/tools -p 'test_*.py' -v
```

- `test_crate_boundaries.py` enforces the eight-crate production workspace,
  one-way normal dependencies, crate-level contract tests, the exact supported
  feature inventory and expansions, workspace MSRV inheritance, explicit
  suppression of downstream default features, and the rule that every Rust
  `cfg(feature = "...")`
  reference is declared by its owning crate. Its lightweight Rust lexer
  ignores comments and string contents while checking that development
  features and retired module names stay outside product source.
- `generate_softclip_truth.py OUTPUT` generates deterministic 151-bp
  directional pairs, known 3' adapter/quality tails, high-quality guaranteed
  5' mismatches, two-ended terminal mismatches, exact origin truth, and
  equal-best decoys for soft-clip ambiguity tests.
- `test_evaluate_mapq_prauc.py` validates tie handling, missing-pair recall
  denominators, operating-point F1, and truth-ledger cardinality checks.
