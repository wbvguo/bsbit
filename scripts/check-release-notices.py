#!/usr/bin/env python3
"""Validate and assemble bsbit's production distribution licenses."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import tomllib
from typing import NoReturn


SCHEMA = "bsbit-third-party-licenses-v2"
ASSEMBLY_SCHEMA = "bsbit-third-party-license-bundle-v1"
ARTIFACTS = ("binary", "source")
BINARY_ARTIFACTS = ("binary",)
APACHE_ASSET = (
    "LICENSE-APACHE",
    "LICENSE-APACHE",
    "a9040321c3712d8fd0b09cf52b17445de04a23a10165049ae187cd39e5c86be5",
)
LICENSE_POLICY_FILES = {
    Path("README.md"),
    Path("THIRD_PARTY_NOTICES.md"),
    Path("license-manifest.json"),
}
AUDITED_REGISTRY = {
    ("block-buffer", "0.12.1"): (
        "MIT OR Apache-2.0",
        "Apache-2.0",
        APACHE_ASSET,
    ),
    ("cfg-if", "1.0.4"): (
        "MIT OR Apache-2.0",
        "Apache-2.0",
        APACHE_ASSET,
    ),
    ("cpufeatures", "0.3.0"): (
        "MIT OR Apache-2.0",
        "Apache-2.0",
        APACHE_ASSET,
    ),
    ("crypto-common", "0.2.2"): (
        "MIT OR Apache-2.0",
        "Apache-2.0",
        APACHE_ASSET,
    ),
    ("digest", "0.11.3"): (
        "MIT OR Apache-2.0",
        "Apache-2.0",
        APACHE_ASSET,
    ),
    ("hybrid-array", "0.4.14"): (
        "MIT OR Apache-2.0",
        "Apache-2.0",
        APACHE_ASSET,
    ),
    ("libc", "0.2.189"): (
        "MIT OR Apache-2.0",
        "Apache-2.0",
        APACHE_ASSET,
    ),
    ("sha2", "0.11.0"): (
        "MIT OR Apache-2.0",
        "Apache-2.0",
        APACHE_ASSET,
    ),
    ("typenum", "1.20.1"): (
        "MIT OR Apache-2.0",
        "Apache-2.0",
        APACHE_ASSET,
    ),
}
AUDITED_NATIVE = {
    "htslib": (
        "1.24",
        "4b705e4fada8ee2b6b15746f725ee8ac51631803",
        "MIT AND BSD-3-Clause",
        "MIT AND BSD-3-Clause",
        (
            "external/htslib/LICENSE",
            "third-party/htslib-1.24-LICENSE",
            "62a9257fa98697f92cf0fee949d59c8afa150b61be210b1f34f5a1bdc6aeb6dd",
        ),
    ),
    "htscodecs": (
        "1.6.7",
        "b9fc194f772e45bb0a1f44b08cbf8697a1384bae",
        "BSD-3-Clause AND LicenseRef-Public-Domain AND CC0-1.0",
        "BSD-3-Clause AND LicenseRef-Public-Domain AND CC0-1.0",
        (
            "external/htslib/htscodecs/LICENSE.md",
            "third-party/htscodecs-1.6.7-LICENSE.md",
            "d17f1bae81abcca5928c2ee0adbc7d4429cc5ad8649d195e011a69f1c3b6b2ad",
        ),
    ),
    "libsais": (
        "2.10.4",
        "ce90878d784b5ff7d019300535675e4a2e22aae0",
        "Apache-2.0",
        "Apache-2.0",
        (
            "external/libsais/LICENSE",
            "third-party/libsais-2.10.4-LICENSE",
            "cfc7749b96f63bd31c3c42b5c471bf756814053e847c10f3eb003417bc523d30",
        ),
    ),
}
AUDITED_RUST_TOOLCHAIN = {
    "rust-standard-library": (
        "1.94.0",
        "4a4ef493e3a1488c6e321570238084b38948f6db",
        "Apache-2.0 OR MIT, plus documented exceptions",
        "Apache-2.0 plus documented exceptions",
        (
            "rust-sysroot/share/doc/rust/COPYRIGHT-library.html",
            "third-party/rust-1.94.0-COPYRIGHT-library.html",
            "af70aaabed1b73e872f14f9130db37e09f3f4d73a5f7c598b9173697a5d2729f",
        ),
    ),
}
AUDITED_SCOPES = {
    ("native", "htslib", "1.24"): ARTIFACTS,
    ("native", "htscodecs", "1.6.7"): ARTIFACTS,
    ("native", "libsais", "2.10.4"): ARTIFACTS,
    (
        "rust-toolchain",
        "rust-standard-library",
        "1.94.0",
    ): BINARY_ARTIFACTS,
    **{
        ("rust-registry", name, version): BINARY_ARTIFACTS
        for name, version in AUDITED_REGISTRY
    },
}


class NoticeError(RuntimeError):
    """A release-license invariant was violated."""


def fail(message: str) -> NoReturn:
    raise NoticeError(message)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json_object(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        fail(f"cannot load {path}: {error}")
    if not isinstance(value, dict):
        fail(f"{path} must contain a JSON object")
    return value


def safe_relative(value: object, field: str) -> Path:
    if not isinstance(value, str) or not value:
        fail(f"{field} must be a non-empty string")
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        fail(f"unsafe {field}: {value!r}")
    return relative


def command_output(arguments: list[str], label: str) -> str:
    try:
        completed = subprocess.run(
            arguments,
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
        )
    except OSError as error:
        fail(f"cannot run {label}: {error}")
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        fail(f"{label} failed with status {completed.returncode}: {detail}")
    return completed.stdout.strip()


def rust_toolchain_metadata() -> dict[str, str]:
    metadata: dict[str, str] = {}
    output = command_output(["rustc", "--version", "--verbose"], "inspect rustc")
    for line in output.splitlines():
        key, separator, value = line.partition(":")
        if separator:
            metadata[key.strip()] = value.strip()
    for required in ("release", "commit-hash", "host"):
        if not metadata.get(required):
            fail(f"rustc metadata lacks {required!r}")
    return metadata


def rust_sysroot() -> Path:
    value = command_output(["rustc", "--print", "sysroot"], "locate rustc sysroot")
    root = Path(value)
    if not root.is_absolute() or not root.is_dir():
        fail(f"rustc returned an invalid sysroot: {value!r}")
    return root


def resolve_policy_source(
    repository: Path, _bundle: Path, relative: Path, kind: str
) -> Path:
    if kind == "rust-toolchain":
        if not relative.parts or relative.parts[0] != "rust-sysroot":
            fail("Rust toolchain license source must begin with rust-sysroot/")
        return rust_sysroot().joinpath(*relative.parts[1:])
    return repository / relative


def registry_lock_packages(repository: Path) -> set[tuple[str, str]]:
    lock_path = repository / "Cargo.lock"
    try:
        lock = tomllib.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, tomllib.TOMLDecodeError) as error:
        fail(f"cannot load {lock_path}: {error}")
    packages = lock.get("package")
    if not isinstance(packages, list):
        fail("Cargo.lock has no package array")

    locked: set[tuple[str, str]] = set()
    for package in packages:
        if not isinstance(package, dict):
            fail("Cargo.lock package entry is not a table")
        source = package.get("source")
        if not isinstance(source, str) or not source.startswith("registry+"):
            continue
        name = package.get("name")
        version = package.get("version")
        if not isinstance(name, str) or not isinstance(version, str):
            fail("registry package lacks a string name/version")
        locked.add((name, version))
    return locked


def registry_source_root(name: str, version: str) -> Path | None:
    cargo_home = Path(os.environ.get("CARGO_HOME", Path.home() / ".cargo"))
    registry_root = cargo_home / "registry" / "src"
    if not registry_root.is_dir():
        return None
    candidates = sorted(registry_root.glob(f"*/{name}-{version}"))
    if len(candidates) > 1:
        fail(f"ambiguous registry source for {name} {version}: {candidates}")
    return candidates[0] if candidates else None


def project_license_ready(repository: Path) -> tuple[bool, str]:
    def nonempty_file(path: Path) -> bool:
        return path.is_file() and path.stat().st_size > 0

    generic_license = any(
        nonempty_file(path)
        for path in (repository / "LICENSE", repository / "LICENSE.md")
    )
    mit_license = nonempty_file(repository / "LICENSE-MIT")
    apache_license = nonempty_file(repository / "LICENSE-APACHE")
    if not generic_license and not mit_license and not apache_license:
        return False, "repository has no project-level license text"

    try:
        root_manifest = tomllib.loads(
            (repository / "Cargo.toml").read_text(encoding="utf-8")
        )
    except (OSError, UnicodeError, tomllib.TOMLDecodeError) as error:
        return False, f"cannot read root Cargo.toml: {error}"
    workspace = root_manifest.get("workspace", {})
    members = workspace.get("members", []) if isinstance(workspace, dict) else []
    workspace_package = (
        workspace.get("package", {}) if isinstance(workspace, dict) else {}
    )
    if not isinstance(members, list) or not isinstance(workspace_package, dict):
        return False, "invalid workspace license metadata"

    def declaration(
        package: dict[str, object], manifest_directory: Path
    ) -> tuple[str, str] | None:
        license_value = package.get("license")
        license_file_value = package.get("license-file")
        if isinstance(license_value, dict) and license_value.get("workspace") is True:
            license_value = workspace_package.get("license")
        if (
            isinstance(license_file_value, dict)
            and license_file_value.get("workspace") is True
        ):
            license_file_value = workspace_package.get("license-file")
            manifest_directory = repository
        if isinstance(license_value, str) and license_value.strip():
            return "expression", license_value.strip()
        if isinstance(license_file_value, str) and license_file_value.strip():
            license_path = manifest_directory / license_file_value
            if nonempty_file(license_path):
                return "file", os.fspath(license_path.resolve())
        return None

    missing: list[str] = []
    declarations: set[tuple[str, str]] = set()
    for member in members:
        if not isinstance(member, str):
            return False, "workspace member is not a string"
        manifest_path = repository / member / "Cargo.toml"
        try:
            manifest = tomllib.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, tomllib.TOMLDecodeError) as error:
            return False, f"cannot read {manifest_path}: {error}"
        package = manifest.get("package", {})
        if not isinstance(package, dict):
            missing.append(member)
            continue
        resolved = declaration(package, repository / member)
        if resolved is None:
            missing.append(str(package.get("name", member)))
        else:
            declarations.add(resolved)
    if missing:
        return False, "Cargo license metadata missing for: " + ", ".join(missing)
    if len(declarations) != 1:
        return False, f"local Cargo license declarations differ: {declarations}"

    kind, value = next(iter(declarations))
    if kind == "file":
        return True, os.fspath(Path(value).relative_to(repository))
    if value == "MIT OR Apache-2.0":
        if not mit_license:
            return False, "MIT OR Apache-2.0 declaration lacks root LICENSE-MIT"
        if not apache_license:
            return False, "MIT OR Apache-2.0 declaration lacks root LICENSE-APACHE"
    elif value == "MIT" and not (generic_license or mit_license):
        return False, "MIT declaration lacks a matching root license file"
    elif value == "Apache-2.0" and not (generic_license or apache_license):
        return False, "Apache-2.0 declaration lacks a matching root license file"
    elif value not in {"MIT", "Apache-2.0"} and not generic_license:
        return False, "Cargo license expression lacks a root LICENSE or LICENSE.md"
    return True, value


def verify_registry_source(
    name: str,
    version: str,
    declared: str,
    selected: str,
    policy_source: Path,
) -> None:
    source_root = registry_source_root(name, version)
    if source_root is None:
        return
    try:
        source_manifest = tomllib.loads(
            (source_root / "Cargo.toml").read_text(encoding="utf-8")
        )
    except (OSError, UnicodeError, tomllib.TOMLDecodeError) as error:
        fail(f"cannot load registry manifest for {name} {version}: {error}")
    source_package = source_manifest.get("package")
    source_license = (
        source_package.get("license") if isinstance(source_package, dict) else None
    )
    if source_license != declared:
        fail(
            f"registry license metadata mismatch for {name} {version}: "
            f"{source_license!r} != {declared!r}"
        )

    notice_files = sorted(
        path.name
        for path in source_root.iterdir()
        if path.is_file() and path.name.upper().startswith("NOTICE")
    )
    if notice_files:
        fail(f"unhandled NOTICE files for {name} {version}: {notice_files}")

    license_files = sorted(
        path
        for path in source_root.iterdir()
        if path.is_file()
        and path.name.upper().startswith(("LICENSE", "COPYING"))
    )
    if selected == "Apache-2.0":
        has_apache = any("APACHE" in path.name.upper() for path in license_files)
        if not has_apache:
            fail(f"registry source lacks Apache license evidence for {name} {version}")
    elif not any(path.read_bytes() == policy_source.read_bytes() for path in license_files):
        fail(
            f"selected {selected} text for {name} {version} does not match "
            "an upstream license file"
        )


def validate(
    repository: Path,
    bundle: Path,
    require_project_license: bool,
    artifact: str | None = None,
) -> list[dict[str, object]]:
    actual_policy_files = {
        path.relative_to(bundle) for path in bundle.rglob("*") if path.is_file()
    }
    if actual_policy_files != LICENSE_POLICY_FILES:
        missing = sorted(
            os.fspath(path) for path in LICENSE_POLICY_FILES - actual_policy_files
        )
        extra = sorted(
            os.fspath(path) for path in actual_policy_files - LICENSE_POLICY_FILES
        )
        fail(f"license policy file set mismatch; missing={missing}, extra={extra}")

    policy = load_json_object(bundle / "license-manifest.json")
    if policy.get("schema") != SCHEMA:
        fail(f"unexpected license schema: {policy.get('schema')!r}")
    components = policy.get("components")
    if not isinstance(components, list) or not components:
        fail("license policy needs a non-empty components array")

    registry_components: dict[tuple[str, str], tuple[str, str, tuple[str, str, str]]] = {}
    native_components: dict[str, tuple[str, str, str, str, tuple[str, str, str]]] = {}
    rust_toolchain_components: dict[
        str, tuple[str, str, str, str, tuple[str, str, str]]
    ] = {}
    observed_scopes: dict[tuple[str, str, str], tuple[str, ...]] = {}
    outputs: dict[Path, tuple[Path, str]] = {}
    for component in components:
        if not isinstance(component, dict):
            fail("component entry is not an object")
        kind = component.get("kind")
        name = component.get("name")
        version = component.get("version")
        declared = component.get("declared_license")
        selected = component.get("selected_license")
        if not all(
            isinstance(value, str)
            for value in (kind, name, version, declared, selected)
        ):
            fail(f"component has non-string identity fields: {component!r}")

        artifact_values = component.get("artifacts")
        if (
            not isinstance(artifact_values, list)
            or not artifact_values
            or not all(isinstance(value, str) for value in artifact_values)
        ):
            fail(f"component {name} has an invalid artifacts array")
        artifacts = tuple(artifact_values)
        if len(set(artifacts)) != len(artifacts):
            fail(f"component {name} repeats an artifact scope")
        unknown_artifacts = sorted(set(artifacts) - set(ARTIFACTS))
        if unknown_artifacts:
            fail(f"component {name} has unknown artifact scopes: {unknown_artifacts}")
        scope_identity = (kind, name, version)
        if scope_identity in observed_scopes:
            fail(f"duplicate scoped component {kind} {name} {version}")
        observed_scopes[scope_identity] = artifacts
        validate_material = artifact is None or artifact in artifacts

        license_file = component.get("license_file")
        if not isinstance(license_file, dict):
            fail(f"component {name} has no license_file object")
        if set(license_file) != {"source", "output", "sha256"}:
            fail(f"component {name} has unexpected license_file fields")
        source_relative = safe_relative(license_file.get("source"), "license source")
        output_relative = safe_relative(license_file.get("output"), "license output")
        expected_hash = license_file.get("sha256")
        if (
            not isinstance(expected_hash, str)
            or len(expected_hash) != 64
            or any(character not in "0123456789abcdef" for character in expected_hash)
        ):
            fail(f"invalid license hash for {name}")
        policy_source: Path | None = None
        if validate_material:
            policy_source = resolve_policy_source(
                repository, bundle, source_relative, kind
            )
            if not policy_source.is_file():
                fail(f"missing authoritative license source {source_relative}")
            observed_hash = file_sha256(policy_source)
            if observed_hash != expected_hash:
                fail(
                    f"license hash mismatch for {source_relative}: "
                    f"{observed_hash} != {expected_hash}"
                )

        asset = (
            source_relative.as_posix(),
            output_relative.as_posix(),
            expected_hash,
        )
        previous = outputs.get(output_relative)
        current = (source_relative, expected_hash)
        if previous is not None and previous != current:
            fail(f"conflicting assembly output {output_relative}")
        outputs[output_relative] = current

        if kind == "rust-registry":
            if set(component) != {
                "kind",
                "name",
                "version",
                "declared_license",
                "selected_license",
                "artifacts",
                "license_file",
            }:
                fail(f"unexpected registry component fields for {name}")
            identity = (name, version)
            if identity in registry_components:
                fail(f"duplicate registry component {name} {version}")
            registry_components[identity] = (declared, selected, asset)
            allowed = {option.strip(" ()") for option in declared.split(" OR ")}
            if selected not in allowed:
                fail(f"selected license {selected!r} is not allowed for {name}")
            if policy_source is not None:
                verify_registry_source(
                    name, version, declared, selected, policy_source
                )
        elif kind == "native":
            if set(component) != {
                "kind",
                "name",
                "version",
                "revision",
                "declared_license",
                "selected_license",
                "artifacts",
                "license_file",
            }:
                fail(f"unexpected native component fields for {name}")
            revision = component.get("revision")
            if not isinstance(revision, str):
                fail(f"native component {name} lacks revision")
            if name in native_components:
                fail(f"duplicate native component {name}")
            native_components[name] = (
                version,
                revision,
                declared,
                selected,
                asset,
            )
        elif kind == "rust-toolchain":
            if set(component) != {
                "kind",
                "name",
                "version",
                "revision",
                "declared_license",
                "selected_license",
                "artifacts",
                "license_file",
            }:
                fail(f"unexpected Rust toolchain component fields for {name}")
            revision = component.get("revision")
            if not isinstance(revision, str):
                fail(f"Rust toolchain component {name} lacks revision")
            if name in rust_toolchain_components:
                fail(f"duplicate Rust toolchain component {name}")
            if validate_material:
                metadata = rust_toolchain_metadata()
                if metadata["release"] != version:
                    fail(
                        f"Rust toolchain release changed: "
                        f"{metadata['release']} != {version}"
                    )
                if metadata["commit-hash"] != revision:
                    fail(
                        f"Rust toolchain revision changed: "
                        f"{metadata['commit-hash']} != {revision}"
                    )
                if metadata["host"] != "x86_64-unknown-linux-gnu":
                    fail(
                        "unsupported audited Rust toolchain host: "
                        f"{metadata['host']}"
                    )
            rust_toolchain_components[name] = (
                version,
                revision,
                declared,
                selected,
                asset,
            )
        else:
            fail(f"unsupported component kind {kind!r}")

    if registry_components != AUDITED_REGISTRY:
        fail("registry license policy differs from audited closure")
    if native_components != AUDITED_NATIVE:
        fail("native license policy differs from audited closure")
    if rust_toolchain_components != AUDITED_RUST_TOOLCHAIN:
        fail("Rust toolchain license policy differs from audited closure")
    if observed_scopes != AUDITED_SCOPES:
        fail("artifact scope policy differs from audited closure")

    locked = registry_lock_packages(repository)
    if locked != set(AUDITED_REGISTRY):
        missing = sorted(set(AUDITED_REGISTRY) - locked)
        extra = sorted(locked - set(AUDITED_REGISTRY))
        fail(
            "Cargo.lock registry closure differs from audited allow-list; "
            f"missing={missing}, extra={extra}"
        )

    notice_text = (bundle / "THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8")
    for required in (
        "Scope matrix",
        "`binary`",
        "`source`",
        "HTSlib",
        "htscodecs",
        "libsais",
        "Rust 1.94.0 standard library",
        "Nine locked Rust registry packages",
        "Apache-2.0",
        "system zlib",
        "Project license",
    ):
        if required not in notice_text:
            fail(f"THIRD_PARTY_NOTICES.md lacks {required!r}")

    ready, reason = project_license_ready(repository)
    if require_project_license and not ready:
        fail(f"project distribution license is not ready: {reason}")
    return components


def components_for_artifact(
    components: list[dict[str, object]], artifact: str
) -> list[dict[str, object]]:
    return [
        component
        for component in components
        if artifact in component["artifacts"]
    ]


def filtered_manifest_bytes(
    artifact: str, components: list[dict[str, object]]
) -> bytes:
    policy = {
        "schema": ASSEMBLY_SCHEMA,
        "artifact": artifact,
        "components": components_for_artifact(components, artifact),
    }
    return (json.dumps(policy, indent=2) + "\n").encode("utf-8")


def artifact_notice_bytes(
    artifact: str, components: list[dict[str, object]]
) -> bytes:
    titles = {
        "binary": "the complete bsbit command-line binary",
        "source": "the recursive bsbit source release",
    }
    selected = components_for_artifact(components, artifact)
    native = [component for component in selected if component["kind"] == "native"]
    registry = [
        component for component in selected if component["kind"] == "rust-registry"
    ]
    toolchain = [
        component for component in selected if component["kind"] == "rust-toolchain"
    ]

    lines = [
        f"# Third-party notices for {titles[artifact]}",
        "",
        "This file was generated from bsbit's audited license manifest. It lists",
        "only components whose bytes are included in this release scope. The",
        "referenced license and copyright files accompany this notice.",
        "",
        "## Native components",
        "",
        "| Component | Version | Selected terms | Included license file |",
        "|---|---|---|---|",
    ]
    for component in native:
        license_file = component["license_file"]
        lines.append(
            f"| {component['name']} | {component['version']} | "
            f"{component['selected_license']} | `{license_file['output']}` |"
        )
    if artifact in BINARY_ARTIFACTS:
        lines.extend(
            [
                "",
                "HTSlib and htscodecs are both retained because linked codec and CRAM",
                "objects are present even though bsbit rejects CRAM input.",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "These components are recursively included in the source release and",
                "retain their complete upstream license texts.",
            ]
        )

    if registry:
        packages = ", ".join(
            f"`{component['name']} {component['version']}`"
            for component in registry
        )
        lines.extend(
            [
                "",
                "## Locked Rust registry packages",
                "",
                packages + ".",
                "",
                "These packages are redistributed under their Apache-2.0 option and",
                "share the single unmodified `LICENSE-APACHE` included at this",
                "directory's root. Their audited sources contain no top-level NOTICE.",
            ]
        )

    if toolchain:
        component = toolchain[0]
        license_file = component["license_file"]
        lines.extend(
            [
                "",
                "## Rust standard library",
                "",
                f"The binary statically includes Rust {component['version']} standard-library",
                "code. Its generated copyright inventory, including applicable",
                "third-party and Unicode exceptions, is included as",
                f"`{license_file['output']}`.",
            ]
        )

    if artifact in BINARY_ARTIFACTS:
        host_libraries = (
            "System zlib, libdeflate, liblzma, libbz2, libgomp, libm, "
            "libgcc_s, libc,"
        )
        lines.extend(
            [
                "",
                "## Dynamically linked host components",
                "",
                host_libraries,
                "and the ELF loader are host prerequisites and are not files in this",
                "standalone binary release. A package or container that ships those libraries",
                "must inventory its exact redistributed versions separately.",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "The source scope includes the recursively pinned native sources. It",
                "does not include Cargo registry source or a compiled Rust standard",
                "library; those materials therefore do not appear in this notice.",
            ]
        )

    lines.extend(
        [
            "",
            "## Project license",
            "",
            "Project-owned work is available under `MIT OR Apache-2.0`; see",
            "`LICENSE-MIT` and `LICENSE-APACHE`. Those terms do not relicense the",
            "components above.",
            "",
        ]
    )
    return "\n".join(lines).encode("utf-8")


def expected_assembly(
    repository: Path,
    bundle: Path,
    components: list[dict[str, object]],
    artifact: str,
) -> dict[Path, bytes]:
    expected = {
        Path("LICENSE-MIT"): (repository / "LICENSE-MIT").read_bytes(),
        Path("LICENSE-APACHE"): (repository / "LICENSE-APACHE").read_bytes(),
        Path("THIRD_PARTY_NOTICES.md"): artifact_notice_bytes(
            artifact, components
        ),
        Path("license-manifest.json"): filtered_manifest_bytes(
            artifact, components
        ),
    }
    for component in components_for_artifact(components, artifact):
        license_file = component["license_file"]
        source_relative = safe_relative(license_file["source"], "license source")
        output_relative = safe_relative(license_file["output"], "license output")
        source = resolve_policy_source(
            repository, bundle, source_relative, component["kind"]
        )
        source_bytes = source.read_bytes()
        previous = expected.get(output_relative)
        if previous is not None and previous != source_bytes:
            fail(f"conflicting bytes for assembly output {output_relative}")
        expected[output_relative] = source_bytes
    return expected


def validate_assembly(
    repository: Path,
    bundle: Path,
    assembly: Path,
    components: list[dict[str, object]],
    artifact: str,
) -> None:
    expected = expected_assembly(repository, bundle, components, artifact)
    actual = {
        path.relative_to(assembly) for path in assembly.rglob("*") if path.is_file()
    }
    if actual != set(expected):
        missing = sorted(os.fspath(path) for path in set(expected) - actual)
        extra = sorted(os.fspath(path) for path in actual - set(expected))
        fail(f"assembled license file set mismatch; missing={missing}, extra={extra}")
    for relative, expected_bytes in expected.items():
        if (assembly / relative).read_bytes() != expected_bytes:
            fail(f"assembled license differs from source: {relative}")


def assemble(
    repository: Path,
    bundle: Path,
    output: Path,
    components: list[dict[str, object]],
    artifact: str,
) -> None:
    output = output.resolve()
    if output.exists():
        fail(f"assembly output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.tmp-", dir=output.parent)
    )
    try:
        for relative, source_bytes in expected_assembly(
            repository, bundle, components, artifact
        ).items():
            destination = temporary / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(source_bytes)
        validate_assembly(repository, bundle, temporary, components, artifact)
        temporary.rename(output)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def expect_failure(
    repository: Path,
    bundle: Path,
    expected_fragment: str,
) -> None:
    try:
        validate(repository, bundle, False)
    except NoticeError as error:
        if expected_fragment not in str(error):
            fail(
                f"negative self-test expected {expected_fragment!r}, got {error!s}"
            )
        return
    fail(f"negative self-test unexpectedly accepted {expected_fragment!r}")


def self_test(repository: Path, bundle: Path) -> None:
    components = validate(repository, bundle, False)
    with tempfile.TemporaryDirectory(prefix="bsbit-license-check-") as temporary:
        temporary_path = Path(temporary)

        mutated = temporary_path / "hash"
        shutil.copytree(bundle, mutated)
        mutated_policy = load_json_object(mutated / "license-manifest.json")
        mutated_components = mutated_policy["components"]
        native = next(
            component
            for component in mutated_components
            if component.get("kind") == "native"
        )
        native["license_file"]["sha256"] = "0" * 64
        (mutated / "license-manifest.json").write_text(
            json.dumps(mutated_policy, indent=2) + "\n", encoding="utf-8"
        )
        expect_failure(repository, mutated, "license hash mismatch")

        extra = temporary_path / "extra"
        shutil.copytree(bundle, extra)
        (extra / "UNDECLARED").write_text("extra\n", encoding="utf-8")
        expect_failure(repository, extra, "license policy file set mismatch")

        drift = temporary_path / "drift"
        shutil.copytree(bundle, drift)
        drift_policy = load_json_object(drift / "license-manifest.json")
        drift_components = drift_policy["components"]
        registry = next(
            component
            for component in drift_components
            if component.get("kind") == "rust-registry"
        )
        registry["version"] = "0.0.0-mutated"
        (drift / "license-manifest.json").write_text(
            json.dumps(drift_policy, indent=2) + "\n", encoding="utf-8"
        )
        expect_failure(repository, drift, "registry license policy")

        scope_drift = temporary_path / "scope-drift"
        shutil.copytree(bundle, scope_drift)
        scope_policy = load_json_object(scope_drift / "license-manifest.json")
        scope_components = scope_policy["components"]
        libsais = next(
            component
            for component in scope_components
            if component.get("name") == "libsais"
        )
        libsais["artifacts"] = ["source"]
        (scope_drift / "license-manifest.json").write_text(
            json.dumps(scope_policy, indent=2) + "\n", encoding="utf-8"
        )
        expect_failure(repository, scope_drift, "artifact scope policy")

        for artifact in ARTIFACTS:
            assembly = temporary_path / f"assembly-{artifact}"
            assemble(repository, bundle, assembly, components, artifact)
            validate_assembly(
                repository, bundle, assembly, components, artifact
            )
            assembled_policy = load_json_object(
                assembly / "license-manifest.json"
            )
            if assembled_policy.get("schema") != ASSEMBLY_SCHEMA:
                fail(f"assembled manifest lost bundle schema for {artifact}")
            if assembled_policy.get("artifact") != artifact:
                fail(f"assembled manifest lost artifact scope {artifact}")
            component_names = {
                component["name"] for component in assembled_policy["components"]
            }
            has_libsais = "libsais" in component_names
            has_rust_std = "rust-standard-library" in component_names
            if not has_libsais:
                fail(f"libsais scope assembly mismatch for {artifact}")
            if has_rust_std != (artifact in BINARY_ARTIFACTS):
                fail(f"Rust standard-library scope mismatch for {artifact}")

        assembly = temporary_path / "assembly-binary"
        (assembly / "third-party/htslib-1.24-LICENSE").write_text(
            "mutation\n", encoding="utf-8"
        )
        try:
            validate_assembly(
                repository, bundle, assembly, components, "binary"
            )
        except NoticeError as error:
            if "assembled license differs" not in str(error):
                fail(f"assembly mutation failed for wrong reason: {error}")
        else:
            fail("assembly mutation was not detected")

        license_repository = temporary_path / "license-repository"
        (license_repository / "crates/a").mkdir(parents=True)
        (license_repository / "crates/b").mkdir(parents=True)
        (license_repository / "Cargo.toml").write_text(
            "[workspace]\n"
            'members = ["crates/a", "crates/b"]\n'
            "[workspace.package]\n"
            'license = "MIT OR Apache-2.0"\n',
            encoding="utf-8",
        )
        inherited_manifest = (
            "[package]\n"
            'name = "{name}"\n'
            'version = "0.0.0"\n'
            "license.workspace = true\n"
        )
        for name in ("a", "b"):
            (license_repository / f"crates/{name}/Cargo.toml").write_text(
                inherited_manifest.format(name=name), encoding="utf-8"
            )
        (license_repository / "LICENSE-MIT").write_text(
            "test MIT text\n", encoding="utf-8"
        )
        apache_path = license_repository / "LICENSE-APACHE"
        apache_path.write_text("test Apache text\n", encoding="utf-8")
        ready, reason = project_license_ready(license_repository)
        if not ready or reason != "MIT OR Apache-2.0":
            fail(f"recommended project-license fixture was rejected: {reason}")
        apache_path.unlink()
        ready, reason = project_license_ready(license_repository)
        if ready or "LICENSE-APACHE" not in reason:
            fail(f"missing Apache license fixture was not rejected: {reason}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle-root", type=Path)
    parser.add_argument("--require-project-license", action="store_true")
    parser.add_argument("--artifact", choices=ARTIFACTS)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--self-test", action="store_true")
    mode.add_argument("--assemble", type=Path, metavar="OUTPUT")
    args = parser.parse_args()
    if args.assemble is not None and args.artifact is None:
        parser.error("--assemble requires --artifact")
    if args.self_test and args.artifact is not None:
        parser.error("--artifact cannot be combined with --self-test")
    return args


def main() -> int:
    args = parse_args()
    repository = Path(__file__).resolve().parent.parent
    bundle = (
        args.bundle_root.resolve()
        if args.bundle_root is not None
        else repository / "external" / "licenses"
    )
    try:
        if args.self_test:
            self_test(repository, bundle)
            components = validate(repository, bundle, args.require_project_license)
        else:
            components = validate(
                repository,
                bundle,
                args.require_project_license or args.assemble is not None,
                args.artifact,
            )
            if args.assemble is not None:
                assemble(
                    repository,
                    bundle,
                    args.assemble,
                    components,
                    args.artifact,
                )
    except (NoticeError, OSError) as error:
        print(f"distribution license validation failed: {error}")
        return 1

    ready, reason = project_license_ready(repository)
    print(
        "third-party license policy OK: "
        f"{len(AUDITED_REGISTRY)} Rust packages, "
        f"{len(AUDITED_NATIVE)} native components, "
        f"{len(AUDITED_RUST_TOOLCHAIN)} Rust toolchain component"
    )
    if args.artifact is not None:
        selected = components_for_artifact(components, args.artifact)
        print(f"artifact scope {args.artifact} OK: {len(selected)} components")
    print(
        f"project distribution license {'OK' if ready else 'PENDING'}: {reason}"
    )
    if args.self_test:
        print(
            "license checker self-test OK: hash, extra-file, policy drift, "
            f"scope drift, {len(ARTIFACTS)} assemblies, project-license readiness"
        )
    if args.assemble is not None:
        print(f"release license directory assembled at {args.assemble.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
