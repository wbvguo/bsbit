"""Repository-level checks for the product crate boundaries."""

from __future__ import annotations

import tomllib
import unittest
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
CRATES = ROOT / "crates"

EXPECTED_CRATES = {
    "bsbit-align",
    "bsbit-call",
    "bsbit-cli",
    "bsbit-combine",
    "bsbit-core",
    "bsbit-hts",
    "bsbit-index",
    "bsbit-io",
}
EXPECTED_RUST_VERSION = "1.89"

# These are dependency ceilings, not requirements. A crate may use fewer edges.
ALLOWED_NORMAL_DEPENDENCIES = {
    "bsbit-core": set(),
    "bsbit-io": set(),
    "bsbit-index": {"bsbit-core", "bsbit-io"},
    "bsbit-hts": {"bsbit-core", "bsbit-io"},
    "bsbit-align": {"bsbit-core", "bsbit-index"},
    "bsbit-call": {"bsbit-core", "bsbit-hts", "bsbit-io"},
    "bsbit-combine": {"bsbit-hts", "bsbit-io"},
    "bsbit-cli": EXPECTED_CRATES - {"bsbit-cli"},
}

DEVELOPMENT_FEATURE_MARKERS = ("experiment", "ablation", "rejected", "superseded")
RETIRED_PRODUCT_PREFIX = "bitmap" + "perbs"
RETIRED_MODULE_NAMES = {f"{RETIRED_PRODUCT_PREFIX}_native", "native", "native_index"}

# Each directory below `src/` must name a real domain rather than an
# implementation status such as "production" or "engine".
EXPECTED_SOURCE_FAMILIES = {
    "bsbit-core": set(),
    "bsbit-io": set(),
    "bsbit-hts": {"bam"},
    "bsbit-index": {"build", "reference", "storage"},
    "bsbit-align": {"paired_end", "search", "single_end", "verification"},
    "bsbit-call": {"evidence", "joint", "meth", "region", "snp"},
    "bsbit-combine": set(),
    "bsbit-cli": {"command", "record_composition"},
}

# These exact implementation boundaries justify one additional directory
# level. No other family may grow deeper without an explicit architecture
# change here.
EXPECTED_NESTED_SOURCE_FAMILIES = {
    "bsbit-align": {Path("verification/narrow")},
    "bsbit-index": {Path("reference/runtime")},
    "bsbit-cli": {Path("command/align")},
}

# Promotion into this inventory is an explicit architecture change. Candidate
# switches are exercised under `agent/worktree/`, not appended opportunistically
# to a product manifest.
EXPECTED_PRODUCT_FEATURE_EXPANSIONS = {
    "bsbit-core": {},
    "bsbit-io": {},
    "bsbit-hts": {},
    "bsbit-call": {},
    "bsbit-combine": {},
    "bsbit-index": {
        "default": [],
        "combined-index": ["dep:libc"],
        "index-construction": ["combined-index"],
    },
    "bsbit-align": {},
    "bsbit-cli": {},
}
EXPECTED_PRODUCT_FEATURES = {
    crate: set(features)
    for crate, features in EXPECTED_PRODUCT_FEATURE_EXPANSIONS.items()
}

# Feature-bearing internal dependencies must not inherit a future default
# capability accidentally. Stable product features opt into every downstream
# capability explicitly through the expansions above.
NO_DEFAULT_FEATURE_DEPENDENCIES = {
    ("bsbit-align", "dependencies", "bsbit-index"),
    ("bsbit-cli", "dependencies", "bsbit-index"),
}

# Keep the index reader and construction closures distinct: alignment may map
# the current index, while only the CLI's index command pulls in its builder.
EXPECTED_INTERNAL_DEPENDENCY_FEATURES = {
    ("bsbit-align", "dependencies", "bsbit-index"): ["combined-index"],
    ("bsbit-cli", "dependencies", "bsbit-index"): ["index-construction"],
}

# These words describe mapping decisions, rather than mechanism-only FM/rank/locate
# capabilities. They must not appear as tokens in the public index API or its
# feature names.
INDEX_POLICY_TOKENS = {
    "alignment",
    "candidate",
    "heuristic",
    "informative",
    "mapq",
    "maximal",
    "mate",
    "pair",
    "policy",
    "rescue",
    "seed",
    "verification",
}

# A feature called "search" would again give index ownership of an aligner
# concern. Exact interval lookup may retain established representation names such as
# `exact_search` and `SearchBase`, so source checks require a policy neighbour
# before treating "search" as a violation.
INDEX_FEATURE_POLICY_TOKENS = INDEX_POLICY_TOKENS | {"search"}
INDEX_PRIVATE_POLICY_TOKENS = {"informative", "maximal", "seed"}

PUBLIC_ITEM = re.compile(
    r"^\s*pub\s+(?:(?:async|const)\s+)?"
    r"(?:fn|struct|enum|trait|type|mod|const)\s+([A-Za-z][A-Za-z0-9_]*)"
)
ANY_ITEM = re.compile(
    r"^\s*(?:pub(?:\([^)]*\))?\s+)?(?:(?:async|const)\s+)?"
    r"(?:fn|struct|enum|trait|type|mod|const)\s+([A-Za-z][A-Za-z0-9_]*)"
)


def identifier_tokens(identifier: str) -> set[str]:
    """Split snake_case and CamelCase identifiers into lowercase words."""

    words: set[str] = set()
    for component in identifier.split("_"):
        words.update(
            word.lower()
            for word in re.findall(
                r"[A-Z]+(?=[A-Z][a-z]|\d|$)|[A-Z]?[a-z]+|\d+", component
            )
        )
    return words


def load_toml(path: Path) -> dict:
    with path.open("rb") as stream:
        return tomllib.load(stream)


def rust_tokens(source: str):
    """Yield code tokens while ignoring comments and string contents.

    String literals are yielded as a single ``string`` token so a real
    ``feature = "..."`` predicate can be recognized without mistaking text in
    comments, documentation, ordinary strings, or raw strings for Rust code.
    The lexer is intentionally small, but handles nested block comments and
    all string prefixes that can contain an apparent ``cfg`` expression.
    """

    index = 0
    length = len(source)
    while index < length:
        if source[index].isspace():
            index += 1
            continue

        if source.startswith("//", index):
            newline = source.find("\n", index + 2)
            index = length if newline == -1 else newline + 1
            continue

        if source.startswith("/*", index):
            depth = 1
            index += 2
            while index < length and depth:
                if source.startswith("/*", index):
                    depth += 1
                    index += 2
                elif source.startswith("*/", index):
                    depth -= 1
                    index += 2
                else:
                    index += 1
            continue

        raw_start = index
        if source.startswith(("br", "cr"), raw_start):
            raw_start += 1
        if raw_start < length and source[raw_start] == "r":
            hash_start = raw_start + 1
            quote = hash_start
            while quote < length and source[quote] == "#":
                quote += 1
            if quote < length and source[quote] == '"':
                hashes = source[hash_start:quote]
                terminator = '"' + hashes
                end = source.find(terminator, quote + 1)
                if end == -1:
                    yield "string", source[quote + 1 :]
                    return
                yield "string", source[quote + 1 : end]
                index = end + len(terminator)
                continue

        quote = index
        if source.startswith(("b\"", "c\""), index):
            quote += 1
        if quote < length and source[quote] == '"':
            value: list[str] = []
            cursor = quote + 1
            while cursor < length:
                if source[cursor] == "\\" and cursor + 1 < length:
                    value.append(source[cursor + 1])
                    cursor += 2
                elif source[cursor] == '"':
                    cursor += 1
                    break
                else:
                    value.append(source[cursor])
                    cursor += 1
            yield "string", "".join(value)
            index = cursor
            continue

        if source[index].isalpha() or source[index] == "_":
            end = index + 1
            while end < length and (source[end].isalnum() or source[end] == "_"):
                end += 1
            yield "identifier", source[index:end]
            index = end
            continue

        yield "punctuation", source[index]
        index += 1


def rust_cfg_features(source: str) -> set[str]:
    """Return feature names used by actual ``cfg``/``cfg_attr`` predicates."""

    tokens = list(rust_tokens(source))
    features: set[str] = set()
    index = 0
    while index < len(tokens):
        kind, value = tokens[index]
        if kind != "identifier" or value not in {"cfg", "cfg_attr"}:
            index += 1
            continue

        cursor = index + 1
        if cursor < len(tokens) and tokens[cursor] == ("punctuation", "!"):
            cursor += 1
        if cursor >= len(tokens) or tokens[cursor] != ("punctuation", "("):
            index += 1
            continue

        depth = 1
        cursor += 1
        while cursor < len(tokens) and depth:
            token_kind, token_value = tokens[cursor]
            if (token_kind, token_value) == ("punctuation", "("):
                depth += 1
            elif (token_kind, token_value) == ("punctuation", ")"):
                depth -= 1
            elif (
                depth > 0
                and (token_kind, token_value) == ("identifier", "feature")
                and cursor + 2 < len(tokens)
                and tokens[cursor + 1] == ("punctuation", "=")
                and tokens[cursor + 2][0] == "string"
            ):
                features.add(tokens[cursor + 2][1])
            cursor += 1
        index = cursor

    return features


class CrateBoundaryTests(unittest.TestCase):
    def test_workspace_contains_only_product_crates(self) -> None:
        workspace = load_toml(ROOT / "Cargo.toml")["workspace"]
        members = {Path(member).name for member in workspace["members"]}
        self.assertEqual(members, EXPECTED_CRATES)

    def test_product_crates_inherit_the_verified_workspace_msrv(self) -> None:
        workspace = load_toml(ROOT / "Cargo.toml")["workspace"]
        self.assertEqual(
            workspace["package"]["rust-version"],
            EXPECTED_RUST_VERSION,
        )
        for crate in EXPECTED_CRATES:
            manifest = load_toml(CRATES / crate / "Cargo.toml")
            with self.subTest(crate=crate):
                self.assertIs(
                    manifest["package"]["rust-version"]["workspace"],
                    True,
                )

    def test_normal_dependencies_follow_the_acyclic_product_graph(self) -> None:
        for crate, allowed in ALLOWED_NORMAL_DEPENDENCIES.items():
            with self.subTest(crate=crate):
                manifest = load_toml(CRATES / crate / "Cargo.toml")
                dependencies = {
                    name
                    for name in manifest.get("dependencies", {})
                    if name.startswith("bsbit-")
                }
                self.assertLessEqual(dependencies, allowed)

    def test_cli_exposes_only_the_umbrella_binary(self) -> None:
        manifest = load_toml(CRATES / "bsbit-cli" / "Cargo.toml")
        self.assertEqual(
            {target["name"] for target in manifest.get("bin", [])},
            {"bsbit"},
        )
        self.assertFalse((CRATES / "bsbit-cli" / "src" / "bin").exists())

    def test_alignment_command_modules_match_the_public_command_names(self) -> None:
        commands = CRATES / "bsbit-cli" / "src" / "command"
        align = commands / "align"
        self.assertTrue((align / "mod.rs").is_file())
        self.assertTrue((align / "single.rs").is_file())
        self.assertTrue((align / "paired.rs").is_file())
        self.assertFalse((commands / "single_end.rs").exists())
        self.assertFalse((commands / "align_general.rs").exists())
        self.assertFalse((commands / "pipeline.rs").exists())
        self.assertFalse((commands / "production_align.rs").exists())
        standard = (align / "mod.rs").read_text(encoding="utf-8")
        dispatch = (commands / "mod.rs").read_text(encoding="utf-8")
        self.assertNotIn("align_general", standard)
        self.assertNotIn("AlignGeneral", dispatch)

    def test_single_end_api_is_not_owned_by_paired_end(self) -> None:
        align = CRATES / "bsbit-align" / "src"
        paired = align / "paired_end"
        single = align / "single_end"
        cli_align = CRATES / "bsbit-cli" / "src" / "command" / "align"
        cli_single = cli_align / "single.rs"
        shared = align / "read_mapping.rs"
        combined_search = align / "search" / "combined_adaptive.rs"

        self.assertFalse((paired / "single.rs").exists())
        self.assertTrue((single / "mapper.rs").is_file())
        self.assertFalse((single / "alignment_pool.rs").exists())
        self.assertFalse((single / "correctness_oracle.rs").exists())
        self.assertTrue((single / "mapq.rs").is_file())
        self.assertTrue((paired / "mapq.rs").is_file())
        self.assertFalse((single / "complete.rs").exists())
        self.assertFalse((align / "alignment_pool.rs").exists())
        self.assertFalse((align / "confidence.rs").exists())
        self.assertFalse((align / "mapq.rs").exists())
        self.assertTrue(shared.is_file())
        self.assertTrue(combined_search.is_file())
        paired_source = "\n".join(
            source.read_text(encoding="utf-8") for source in paired.rglob("*.rs")
        )
        for identifier in (
            "SingleAlignmentResult",
            "SingleBatchAligner",
            "SingleMappingStatus",
        ):
            with self.subTest(identifier=identifier):
                self.assertNotIn(identifier, paired_source)

        cli_source = cli_single.read_text(encoding="utf-8")
        self.assertIn("bsbit_align::single_end", cli_source)
        self.assertNotIn("bsbit_align::paired_end", cli_source)
        self.assertNotIn("super::index", cli_source)

        paired_cli_source = "\n".join(
            source.read_text(encoding="utf-8") for source in cli_align.glob("*.rs")
        )
        self.assertNotIn("super::index", paired_cli_source)

        single_source = (single / "mapper.rs").read_text(encoding="utf-8")
        shared_source = shared.read_text(encoding="utf-8")
        combined_search_source = combined_search.read_text(encoding="utf-8")
        for source_name, source in (
            ("single-end mapper", single_source),
            ("shared read mapping", shared_source),
            ("combined search", combined_search_source),
        ):
            with self.subTest(source=source_name):
                self.assertNotIn("crate::paired_end", source)

        self.assertIn("CombinedTwoLaneSearchState", combined_search_source)
        self.assertNotIn("CombinedAdaptiveSearchState", combined_search_source)
        self.assertNotIn("candidate_pair", combined_search_source)
        for declaration in (
            "struct ReadCandidate",
            "struct ReadWorkspace",
            "struct ReadAlignmentMetrics",
        ):
            with self.subTest(declaration=declaration):
                self.assertIn(declaration, shared_source)
                self.assertNotIn(declaration, paired_source)

        placement = (align / "placement.rs").read_text(encoding="utf-8")
        materialize = (align / "materialize.rs").read_text(encoding="utf-8")
        self.assertIn("struct ReadPlacement", placement)
        self.assertIn("fn traceback_read_placement", materialize)
        self.assertNotIn("fn traceback_read_placement", paired_source)

    def test_core_domain_values_are_not_reexported_by_downstream_crates(self) -> None:
        for crate in EXPECTED_CRATES - {"bsbit-core"}:
            for source in (CRATES / crate / "src").rglob("*.rs"):
                with self.subTest(source=source.relative_to(ROOT)):
                    self.assertNotRegex(
                        source.read_text(encoding="utf-8"),
                        r"(?m)^\s*pub\s+use\s+bsbit_core\b",
                        "shared domain values must be imported from bsbit-core directly",
                    )

    def test_every_product_crate_has_contract_tests(self) -> None:
        for crate in EXPECTED_CRATES:
            with self.subTest(crate=crate):
                tests = CRATES / crate / "tests"
                self.assertTrue(tests.is_dir(), f"{crate} has no tests/ directory")
                self.assertTrue(
                    any(tests.rglob("*.rs")),
                    f"{crate} has no crate-level contract test",
                )

    def test_source_trees_use_only_intentional_family_levels(self) -> None:
        for crate, expected_families in EXPECTED_SOURCE_FAMILIES.items():
            src = CRATES / crate / "src"
            observed_families = {
                child.name for child in src.iterdir() if child.is_dir()
            }
            with self.subTest(crate=crate):
                self.assertEqual(observed_families, expected_families)
                for source in src.rglob("*.rs"):
                    relative = source.relative_to(src)
                    if len(relative.parts) <= 2:
                        continue
                    self.assertEqual(len(relative.parts), 3)
                    self.assertIn(
                        relative.parent,
                        EXPECTED_NESTED_SOURCE_FAMILIES.get(crate, set()),
                        f"source tree uses an unapproved nested family: {relative}",
                    )

    def test_development_features_do_not_enter_product_manifests(self) -> None:
        for crate in EXPECTED_CRATES:
            manifest = load_toml(CRATES / crate / "Cargo.toml")
            expansions = manifest.get("features", {})
            features = set(expansions)
            with self.subTest(crate=crate, inventory="promoted"):
                self.assertEqual(features, EXPECTED_PRODUCT_FEATURES[crate])
            with self.subTest(crate=crate, inventory="exact-expansions"):
                self.assertEqual(
                    expansions,
                    EXPECTED_PRODUCT_FEATURE_EXPANSIONS[crate],
                )
            for feature in features:
                with self.subTest(crate=crate, feature=feature):
                    self.assertFalse(
                        any(marker in feature for marker in DEVELOPMENT_FEATURE_MARKERS),
                        f"development feature {feature!r} belongs under agent/worktree",
                    )

    def test_feature_bearing_internal_dependencies_disable_defaults(self) -> None:
        for crate, section, dependency in NO_DEFAULT_FEATURE_DEPENDENCIES:
            manifest = load_toml(CRATES / crate / "Cargo.toml")
            declaration = manifest[section][dependency]
            with self.subTest(
                crate=crate,
                section=section,
                dependency=dependency,
            ):
                self.assertIsInstance(declaration, dict)
                self.assertIs(declaration.get("default-features"), False)

    def test_internal_dependency_features_preserve_reader_builder_separation(self) -> None:
        for key, expected in EXPECTED_INTERNAL_DEPENDENCY_FEATURES.items():
            crate, section, dependency = key
            manifest = load_toml(CRATES / crate / "Cargo.toml")
            declaration = manifest[section][dependency]
            with self.subTest(
                crate=crate,
                section=section,
                dependency=dependency,
            ):
                self.assertEqual(declaration.get("features", []), expected)

    def test_every_cfg_feature_is_declared_by_its_crate(self) -> None:
        for crate in EXPECTED_CRATES:
            crate_root = CRATES / crate
            manifest = load_toml(crate_root / "Cargo.toml")
            declared = set(manifest.get("features", {}))
            references: dict[str, list[Path]] = {}
            for source in crate_root.rglob("*.rs"):
                for feature in rust_cfg_features(source.read_text(encoding="utf-8")):
                    references.setdefault(feature, []).append(source.relative_to(ROOT))

            undeclared = {
                feature: paths
                for feature, paths in references.items()
                if feature not in declared
            }
            with self.subTest(crate=crate):
                details = "; ".join(
                    f"{feature}: {', '.join(map(str, paths))}"
                    for feature, paths in sorted(undeclared.items())
                )
                self.assertFalse(
                    undeclared,
                    f"cfg features missing from {crate}/Cargo.toml: {details}",
                )

    def test_cfg_feature_scanner_ignores_comments_and_string_contents(self) -> None:
        source = r'''
            #[cfg(all(feature = "real-one", not(feature="real-two")))]
            const ENABLED: bool = cfg!(feature = r#"real-three"#);
            #[cfg_attr(feature = "real-four", allow(dead_code))]
            const TEXT: &str = r#"cfg(feature = "string-only")"#;
            // #[cfg(feature = "line-comment")]
            /* #[cfg(feature = "block-comment")]
               /* #[cfg(feature = "nested-comment")] */ */
        '''
        self.assertEqual(
            rust_cfg_features(source),
            {"real-one", "real-two", "real-three", "real-four"},
        )

    def test_retired_experiment_names_do_not_reenter_src(self) -> None:
        for crate in EXPECTED_CRATES:
            src = CRATES / crate / "src"
            self.assertFalse((CRATES / crate / "experiments").exists())
            for source in src.rglob("*.rs"):
                with self.subTest(source=source.relative_to(ROOT)):
                    self.assertNotIn(source.stem.lower(), RETIRED_MODULE_NAMES)
                    self.assertNotIn(
                        RETIRED_PRODUCT_PREFIX,
                        source.read_text(encoding="utf-8").lower(),
                    )

    def test_retired_alignment_parallel_surfaces_do_not_reenter_src(self) -> None:
        retired = {
            "CombinedSearchHeuristicStageProfile",
            "DiscoveryStageProfile",
            "FusedShadow",
            "MateFirstJointRoute",
            "CertifiedPairRoute",
            "ProfiledMate",
            "PairedLibraryProfile::Pbat",
            "controlled A/B",
            "audit-only",
            "profile_stages",
        }
        for source in (CRATES / "bsbit-align" / "src").rglob("*.rs"):
            contents = source.read_text(encoding="utf-8")
            with self.subTest(source=source.relative_to(ROOT)):
                for marker in retired:
                    self.assertNotIn(
                        marker,
                        contents,
                        "unselected alignment candidates belong in agent/worktree or tests",
                    )

    def test_cli_and_index_do_not_expose_oracle_switches(self) -> None:
        cli_forbidden = {
            "--debug",
            "--reference-backend",
            "--global-index-builder",
            "--builder",
            "--audit-samples",
            "--sa-stride",
            "direct-libsais",
            "direct-libsais64",
        }
        for source in (CRATES / "bsbit-cli" / "src").rglob("*.rs"):
            if source.name == "tests.rs":
                continue
            contents = source.read_text(encoding="utf-8")
            with self.subTest(source=source.relative_to(ROOT)):
                for marker in cli_forbidden:
                    self.assertNotIn(
                        marker,
                        contents,
                        "production CLI must not expose debug/oracle selection",
                    )

        index_forbidden = {"BSBIT_INTERNAL_", "BSBIT_REFERENCE_LOAD_PROFILE"}
        for source in (CRATES / "bsbit-index" / "src").rglob("*.rs"):
            contents = source.read_text(encoding="utf-8")
            with self.subTest(source=source.relative_to(ROOT)):
                for marker in index_forbidden:
                    self.assertNotIn(
                        marker,
                        contents,
                        "production index must not have environment-selected behavior",
                    )

        storage_forbidden = {
            "_profiled",
            "advise_reference_huge_pages",
            "collapse_reference_huge_pages",
        }
        for source in (CRATES / "bsbit-index" / "src" / "storage").rglob("*.rs"):
            contents = source.read_text(encoding="utf-8")
            with self.subTest(source=source.relative_to(ROOT)):
                for marker in storage_forbidden:
                    self.assertNotIn(
                        marker,
                        contents,
                        "runtime storage must not expose profiling or huge-page tuning",
                    )

        retired_public_types = {
            "CombinedIndexBuildBackend",
            "GlobalPackedCacheBuilder",
            "Libsais64GlobalBwtBuilder",
            "Libsais32GapGlobalBwtBuilder",
            "Libsais32GapCompactGlobalBwtBuilder",
            "IndependentGlobalBwtProcessBuilder",
        }
        for source in (CRATES / "bsbit-index" / "src").rglob("*.rs"):
            contents = source.read_text(encoding="utf-8")
            with self.subTest(source=source.relative_to(ROOT)):
                for type_name in retired_public_types:
                    self.assertNotRegex(
                        contents,
                        rf"(?m)^pub\\s+(?:struct|enum)\\s+{type_name}\\b",
                        "production index must not export oracle builder selection",
                    )

    def test_index_public_api_and_features_do_not_own_mapping_policy(self) -> None:
        manifest = load_toml(CRATES / "bsbit-index" / "Cargo.toml")
        for feature in manifest.get("features", {}):
            with self.subTest(feature=feature):
                self.assertTrue(
                    INDEX_FEATURE_POLICY_TOKENS.isdisjoint(feature.split("-")),
                    f"mapping-policy feature {feature!r} belongs in bsbit-align",
                )

        for source in (CRATES / "bsbit-index" / "src").rglob("*.rs"):
            for line_number, line in enumerate(
                source.read_text(encoding="utf-8").splitlines(), start=1
            ):
                match = PUBLIC_ITEM.match(line)
                if match is None:
                    continue
                identifier = match.group(1)
                with self.subTest(
                    source=source.relative_to(ROOT),
                    line=line_number,
                    identifier=identifier,
                ):
                    self.assertTrue(
                        INDEX_POLICY_TOKENS.isdisjoint(identifier_tokens(identifier)),
                        f"mapping-policy API {identifier!r} belongs in bsbit-align",
                    )

    def test_index_private_implementation_does_not_reintroduce_seed_policy(self) -> None:
        for source in (CRATES / "bsbit-index" / "src").rglob("*.rs"):
            for line_number, line in enumerate(
                source.read_text(encoding="utf-8").splitlines(), start=1
            ):
                match = ANY_ITEM.match(line)
                if match is None:
                    continue
                identifier = match.group(1)
                tokens = identifier_tokens(identifier)
                with self.subTest(
                    source=source.relative_to(ROOT),
                    line=line_number,
                    identifier=identifier,
                ):
                    self.assertTrue(
                        INDEX_PRIVATE_POLICY_TOKENS.isdisjoint(tokens),
                        f"index implementation policy {identifier!r} belongs in bsbit-align",
                    )
                    self.assertFalse(
                        "search" in tokens
                        and not INDEX_POLICY_TOKENS.isdisjoint(tokens),
                        f"index search policy {identifier!r} belongs in bsbit-align",
                    )

    def test_hts_sam_and_bam_share_the_alignment_record_model(self) -> None:
        hts = CRATES / "bsbit-hts" / "src"
        self.assertTrue((hts / "alignment_record.rs").is_file())
        alignment = (hts / "alignment_record.rs").read_text(encoding="utf-8")
        for declaration in (
            "enum AlignmentCigarOp",
            "struct AlignmentCigarRun",
            "struct BorrowedAlignmentRead",
            "struct BorrowedAlignmentRecord",
            "struct AlignmentRecordBatch",
            "struct AlignmentPlacement",
        ):
            with self.subTest(declaration=declaration):
                self.assertIn(declaration, alignment)
        self.assertNotIn("crate::bam", alignment)

        for codec in ("sam.rs", "bam/writer.rs"):
            with self.subTest(codec=codec):
                source = (hts / codec).read_text(encoding="utf-8")
                self.assertIn("crate::alignment_record", source)

        sam = (hts / "sam.rs").read_text(encoding="utf-8")
        bam = "\n".join(
            source.read_text(encoding="utf-8")
            for source in (hts / "bam").glob("*.rs")
        )
        self.assertIn("BorrowedAlignmentRecord", sam)
        self.assertIn("sam_borrowed_record_bytes", sam)
        self.assertIn("BorrowedAlignmentRecord", bam)
        self.assertIn("write_borrowed_alignment_record", bam)
        for stale_model in (
            "struct " + "Direct" + "AlignmentRecord",
            "struct " + "Direct" + "AlignmentBatch",
            "struct " + "Slab" + "AlignmentRead",
            "struct " + "Ungapped" + "Placement",
        ):
            with self.subTest(stale_model=stale_model):
                self.assertNotIn(stale_model, bam)

        io_source_names = {path.name for path in (CRATES / "bsbit-io" / "src").glob("*.rs")}
        self.assertNotIn("alignment.rs", io_source_names)
        self.assertNotIn("sam.rs", io_source_names)
        self.assertNotIn("bam.rs", io_source_names)

    def test_retired_compact_alignment_api_does_not_reenter_product_crates(self) -> None:
        retired = (
            "Direct" + "AlignmentRecord",
            "Direct" + "AlignmentBatch",
            "Slab" + "AlignmentRead",
            "Ungapped" + "Placement",
            "Direct" + "BamAuxiliaryMode",
            "Direct" + "CigarOp",
            "Direct" + "CigarRun",
            "Narrow" + "CompactPlacementDistances",
            "write_borrowed_" + "direct_bam",
            "write_record_" + "direct_bam",
        )
        for source in CRATES.rglob("*.rs"):
            contents = source.read_text(encoding="utf-8")
            for identifier in retired:
                with self.subTest(source=source.relative_to(ROOT), identifier=identifier):
                    self.assertNotIn(identifier, contents)

    def test_removed_internal_scaffolding_does_not_reenter_product_crates(self) -> None:
        retired = (
            "NormalizedSequence" + "Backing",
            "LoadedReferenceCatalog" + "Snapshot",
            "CombinedIndexBuild" + "Summary",
            "CombinedIndexBuild" + "Timings",
            "load_reference_catalog_" + "snapshot",
            "visit_raw_intervals_wavefront_" + "complete",
        )
        for source in CRATES.rglob("*.rs"):
            contents = source.read_text(encoding="utf-8")
            for identifier in retired:
                with self.subTest(source=source.relative_to(ROOT), identifier=identifier):
                    self.assertNotIn(identifier, contents)

    def test_retired_reference_snapshot_format_stays_out_of_product_source(self) -> None:
        storage = CRATES / "bsbit-index" / "src" / "storage"
        self.assertFalse((storage / "snapshot.rs").exists())
        self.assertFalse((ROOT / "docs" / "reference-snapshot-format.md").exists())
        for source in CRATES.rglob("*.rs"):
            contents = source.read_text(encoding="utf-8")
            with self.subTest(source=source.relative_to(ROOT)):
                self.assertNotIn("ReferenceSnapshot", contents)
                self.assertNotIn("BSBITSNP", contents)

    def test_canonical_alignment_cannot_construct_an_index(self) -> None:
        source = "\n".join(
            path.read_text(encoding="utf-8")
            for path in (CRATES / "bsbit-cli" / "src" / "command" / "align").glob("*.rs")
        )
        self.assertNotIn("bsbit_index::build", source)
        self.assertNotIn("build_combined_index", source)
        self.assertNotIn("load_or_build", source)

    def test_call_modes_and_shared_evidence_have_explicit_ownership(self) -> None:
        call = CRATES / "bsbit-call" / "src"
        self.assertFalse((call / "calling").exists())
        for family in ("meth", "snp", "joint", "evidence", "region"):
            with self.subTest(family=family):
                self.assertTrue((call / family / "mod.rs").is_file())
        contract = (call / "region" / "mod.rs").read_text(encoding="utf-8")
        planner = call / "region" / "planner.rs"
        self.assertTrue(planner.is_file())
        self.assertNotIn("crate::calling", contract)
        self.assertNotIn("struct BamReference", contract)
        self.assertNotIn("struct CallRegion", contract)
        planner_source = planner.read_text(encoding="utf-8")
        self.assertIn("fn plan_call_regions", planner_source)

    def test_combine_input_and_output_contracts_are_not_owned_by_merge(self) -> None:
        combine = CRATES / "bsbit-combine" / "src"
        for name in ("request.rs", "result.rs", "site.rs"):
            with self.subTest(file=name):
                self.assertTrue((combine / name).is_file())
        for name in ("input.rs", "output.rs"):
            source = (combine / name).read_text(encoding="utf-8")
            with self.subTest(file=name):
                self.assertNotIn("crate::merge", source)


if __name__ == "__main__":
    unittest.main()
