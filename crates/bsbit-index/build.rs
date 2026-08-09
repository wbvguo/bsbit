//! Reproducible private build of the pinned libsais source for index builders.

use std::env;
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn main() {
    println!("cargo:rerun-if-env-changed=CARGO_FEATURE_INDEX_CONSTRUCTION");
    if env::var_os("CARGO_FEATURE_INDEX_CONSTRUCTION").is_none() {
        return;
    }

    assert!(
        env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("linux"),
        "native suffix-array builders currently support Linux only"
    );
    assert!(
        env::var("CARGO_CFG_TARGET_ENDIAN").as_deref() == Ok("little"),
        "the combined-index builder requires a little-endian target"
    );

    let manifest_dir = PathBuf::from(required_env("CARGO_MANIFEST_DIR"));
    let repository_root = manifest_dir
        .parent()
        .and_then(Path::parent)
        .expect("index-builder crate must remain below repository/crates")
        .to_path_buf();
    let source_root = repository_root.join("external/libsais");
    assert!(
        source_root.join("include/libsais.h").is_file()
            && source_root.join("src/libsais.c").is_file()
            && source_root.join("src/libsais64.c").is_file(),
        "libsais source is missing; run `git submodule update --init --recursive`"
    );

    println!("cargo:rerun-if-changed={}", source_root.display());
    println!("cargo:rerun-if-env-changed=CC");
    println!("cargo:rerun-if-env-changed=AR");

    let output_root = PathBuf::from(required_env("OUT_DIR"));
    let compiler = env::var("CC").unwrap_or_else(|_| String::from("cc"));
    let archiver = env::var("AR").unwrap_or_else(|_| String::from("ar"));
    let include = source_root.join("include");
    let sources = [
        source_root.join("src/libsais.c"),
        source_root.join("src/libsais64.c"),
    ];
    let mut objects = Vec::new();
    for source in sources {
        let stem = source
            .file_stem()
            .and_then(OsStr::to_str)
            .expect("libsais source stem is UTF-8");
        let object = output_root.join(format!("{stem}.o"));
        run(
            Command::new(&compiler)
                .args([
                    "-std=c11",
                    "-O3",
                    "-DNDEBUG",
                    "-fPIC",
                    "-fopenmp",
                    "-DLIBSAIS_OPENMP",
                ])
                .arg("-I")
                .arg(&include)
                .arg("-c")
                .arg(&source)
                .arg("-o")
                .arg(&object),
            "compile libsais",
        );
        objects.push(object);
    }

    let archive = output_root.join("libbsbit_libsais.a");
    let mut archive_command = Command::new(archiver);
    archive_command.arg("crs").arg(&archive).args(&objects);
    run(&mut archive_command, "archive libsais");
    println!("cargo:rustc-link-search=native={}", output_root.display());
    println!("cargo:rustc-link-lib=static=bsbit_libsais");
    println!("cargo:rustc-link-lib=dylib=gomp");
}

fn required_env(name: &str) -> String {
    env::var(name).unwrap_or_else(|_| panic!("required Cargo variable {name} is missing"))
}

fn run(command: &mut Command, label: &str) {
    let output = run_output(command, label);
    assert!(
        output.status.success(),
        "{label} failed with {}\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn run_output(command: &mut Command, label: &str) -> Output {
    command
        .output()
        .unwrap_or_else(|error| panic!("{label}: {error}"))
}
