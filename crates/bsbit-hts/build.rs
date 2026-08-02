//! Reproducible private build of the pinned `HTSlib` source and project C shim.

use std::env;
use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const NATIVE_SANITIZER_ENV: &str = "BSBIT_NATIVE_SANITIZER";
const THREAD_SANITIZER_CFLAGS: &[&str] = &[
    "-O1",
    "-g",
    "-fno-omit-frame-pointer",
    "-fPIC",
    "-fsanitize=thread",
];

fn main() {
    assert!(
        env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("linux"),
        "bsbit-hts currently supports only the audited Linux build profile"
    );

    let manifest_dir = PathBuf::from(required_env("CARGO_MANIFEST_DIR"));
    let repository_root = manifest_dir
        .parent()
        .and_then(Path::parent)
        .expect("HTS crate must remain below repository/crates")
        .to_path_buf();
    let external = repository_root.join("external/htslib");
    let shim = manifest_dir.join("htslib-shim");
    assert!(
        external.join("hts.c").is_file()
            && external.join("htscodecs/htscodecs/htscodecs.c").is_file(),
        "HTSlib source is missing; run `git submodule update --init --recursive`"
    );

    println!("cargo:rerun-if-changed={}", external.display());
    println!("cargo:rerun-if-changed={}", shim.display());
    println!("cargo:rerun-if-env-changed=CC");
    println!("cargo:rerun-if-env-changed=AR");
    println!("cargo:rerun-if-env-changed=NUM_JOBS");
    println!("cargo:rerun-if-env-changed={NATIVE_SANITIZER_ENV}");

    let native_sanitizer = NativeSanitizer::from_environment();
    let compiler = native_sanitizer.compiler();
    let archiver = env::var("AR").unwrap_or_else(|_| String::from("ar"));

    let output_root = PathBuf::from(required_env("OUT_DIR")).join("native");
    let staged_source = output_root.join("htslib-1.24");
    if output_root.exists() {
        fs::remove_dir_all(&output_root).expect("remove previous private native build");
    }
    copy_tree_without_git(&external, &staged_source);
    let htscodecs_version = staged_source.join("htscodecs/htscodecs/version.h");
    fs::write(
        htscodecs_version,
        "#define HTSCODECS_VERSION_TEXT \"1.6.7\"\n",
    )
    .expect("write audited staged htscodecs version header");

    run(
        Command::new("autoreconf")
            .arg("-i")
            .current_dir(&staged_source),
        "autoreconf HTSlib",
    );
    let mut configure = Command::new("./configure");
    configure
        .arg("--disable-plugins")
        .arg("--disable-libcurl")
        .arg("--disable-gcs")
        .arg("--disable-s3")
        // Keep the native dependency ISA contract aligned with the Rust
        // x86-64-v3 build: retain SSE4/AVX2, but do not compile HTScodecs'
        // independently dispatched AVX-512 CRAM implementation.
        .env(
            "hts_cv_check_cflags_needed_avx512f___mavx512f__mpopcnt",
            "unsupported",
        )
        .current_dir(&staged_source);
    native_sanitizer.configure(&mut configure, &compiler);
    run(&mut configure, "configure HTSlib");
    let jobs = env::var("NUM_JOBS").unwrap_or_else(|_| String::from("1"));
    run(
        Command::new("make")
            .arg(format!("-j{jobs}"))
            .arg("lib-static")
            .current_dir(&staged_source),
        "build static HTSlib",
    );

    let shim_object = output_root.join("bsbit_hts.o");
    let shim_archive = output_root.join("libbsbit_htslib_shim.a");
    let mut compile_shim = Command::new(&compiler);
    compile_shim.args([
        "-std=c11",
        "-D_XOPEN_SOURCE=700",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        "-Wconversion",
        "-Werror",
    ]);
    native_sanitizer.compile(&mut compile_shim);
    compile_shim
        .arg(format!("-I{}", shim.display()))
        .arg("-isystem")
        .arg(&staged_source)
        .arg("-c")
        .arg(shim.join("bsbit_hts.c"))
        .arg("-o")
        .arg(&shim_object);
    run(&mut compile_shim, "compile project HTSlib shim");
    run(
        Command::new(archiver)
            .arg("crs")
            .arg(&shim_archive)
            .arg(&shim_object),
        "archive project HTSlib shim",
    );

    emit_link_directives(&output_root, &staged_source);
}

fn emit_link_directives(output_root: &Path, staged_source: &Path) {
    println!("cargo:rustc-link-search=native={}", output_root.display());
    println!("cargo:rustc-link-search=native={}", staged_source.display());
    println!("cargo:rustc-link-lib=static=bsbit_htslib_shim");
    println!("cargo:rustc-link-lib=static=hts");
    for library in ["deflate", "lzma", "bz2", "z", "m", "pthread"] {
        println!("cargo:rustc-link-lib={library}");
    }
}

#[derive(Clone, Copy)]
enum NativeSanitizer {
    None,
    Thread,
}

impl NativeSanitizer {
    fn from_environment() -> Self {
        match env::var(NATIVE_SANITIZER_ENV) {
            Err(env::VarError::NotPresent) => Self::None,
            Ok(value) if value == "thread" => Self::Thread,
            Ok(value) => panic!(
                "unsupported {NATIVE_SANITIZER_ENV} value {value:?}; expected thread or unset"
            ),
            Err(env::VarError::NotUnicode(_)) => {
                panic!("{NATIVE_SANITIZER_ENV} must be valid Unicode")
            }
        }
    }

    fn compiler(self) -> String {
        match self {
            Self::None => env::var("CC").unwrap_or_else(|_| String::from("cc")),
            Self::Thread => {
                let compiler = env::var("CC").unwrap_or_else(|_| String::from("clang"));
                let executable = Path::new(&compiler)
                    .file_name()
                    .and_then(OsStr::to_str)
                    .unwrap_or("");
                assert!(
                    executable == "clang" || executable.starts_with("clang-"),
                    "{NATIVE_SANITIZER_ENV}=thread requires Clang; observed CC={compiler:?}"
                );
                compiler
            }
        }
    }

    fn configure(self, command: &mut Command, compiler: &str) {
        if matches!(self, Self::Thread) {
            command
                .env("CC", compiler)
                .env("CFLAGS", THREAD_SANITIZER_CFLAGS.join(" "))
                .env("LDFLAGS", "-fsanitize=thread");
            println!("cargo:warning=native ThreadSanitizer instrumentation enabled");
        }
    }

    fn compile(self, command: &mut Command) {
        match self {
            Self::None => {
                command.args(["-O2", "-fPIC"]);
            }
            Self::Thread => {
                command.args(THREAD_SANITIZER_CFLAGS);
            }
        }
    }
}

fn required_env(name: &str) -> String {
    env::var(name).unwrap_or_else(|_| panic!("required Cargo variable {name} is missing"))
}

fn copy_tree_without_git(source: &Path, destination: &Path) {
    fs::create_dir_all(destination).expect("create private native source directory");
    for entry in fs::read_dir(source).expect("read pinned native source directory") {
        let entry = entry.expect("read pinned native source entry");
        if entry.file_name() == OsStr::new(".git") {
            continue;
        }
        let source_path = entry.path();
        let destination_path = destination.join(entry.file_name());
        let metadata = fs::symlink_metadata(&source_path).expect("inspect native source entry");
        if metadata.is_dir() {
            copy_tree_without_git(&source_path, &destination_path);
        } else if metadata.file_type().is_symlink() {
            copy_symlink(&source_path, &destination_path);
        } else {
            fs::copy(&source_path, &destination_path).expect("copy native source file");
        }
    }
}

#[cfg(unix)]
fn copy_symlink(source: &Path, destination: &Path) {
    let target = fs::read_link(source).expect("read native source symlink");
    std::os::unix::fs::symlink(target, destination).expect("copy native source symlink");
}

#[cfg(not(unix))]
fn copy_symlink(_source: &Path, _destination: &Path) {
    panic!("the audited native source copy requires Unix symlink support");
}

fn run(command: &mut Command, label: &str) {
    let output = run_output(command, label);
    if !output.stdout.is_empty() {
        println!("cargo:warning={label} completed");
    }
}

fn run_output(command: &mut Command, label: &str) -> Output {
    let output = command
        .output()
        .unwrap_or_else(|error| panic!("could not execute {label}: {error}"));
    assert!(
        output.status.success(),
        "{label} failed with {}\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}
