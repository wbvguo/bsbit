# Production distribution licenses

This directory contains the audited policy used to assemble license material
for a concrete bsbit release artifact. It sits beside the pinned native
sources because it describes external production dependencies; product code
does not read it at runtime.

The policy has two explicit scopes:

| Scope | Distributed content |
|---|---|
| `binary` | The single `bsbit` binary, including every public subcommand |
| `source` | Project source plus the recursively included native sources |

- `THIRD_PARTY_NOTICES.md` explains the complete scope matrix.
- `license-manifest.json` binds every component to its applicable scopes,
  selected license text, version/revision, and SHA-256 digest.
- Native license texts remain authoritative beside their pinned sources.
- Locked Rust registry packages reuse the unmodified root `LICENSE-APACHE`.
- The binary scope copies the copyright inventory from the exact audited Rust
  standard library installed with Rust 1.94.0. A different toolchain fails
  validation until its identity and generated inventory have been reviewed.

Validate the complete policy and its negative tests with:

```sh
git submodule update --init --recursive
python3 scripts/check-release-notices.py --require-project-license
python3 scripts/check-release-notices.py --self-test
```

Assemble a minimal distributable directory under ignored `dist/` by naming the
artifact explicitly:

```sh
python3 scripts/check-release-notices.py \
  --artifact binary --assemble dist/bsbit/licenses
python3 scripts/check-release-notices.py \
  --artifact source --assemble dist/bsbit-source/licenses
```

The output path must not already exist. Each output receives a filtered
manifest and generated notice containing only components distributed in that
scope. Release archives should include the assembled directory beside their
the `bsbit` binary or source tree.

System libraries that remain dynamically linked are platform prerequisites,
not bytes in the standalone artifact. A package or container that also ships
those libraries must extend the inventory for the exact files and versions it
redistributes.

Experimental dependencies belong to dated records under `agent/worktree/` and
are deliberately absent from every production scope.
