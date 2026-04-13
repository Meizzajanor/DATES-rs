use std::fs;
use std::path::{Path, PathBuf};

use assert_cmd::Command;

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/toy")
}

fn golden_dir() -> PathBuf {
    fixture_dir().join("golden")
}

fn copy_files(names: &[&str], destination: &Path) {
    let fixture = fixture_dir();
    let golden = golden_dir();
    for name in names {
        let source = fixture.join(name);
        let source = if source.exists() {
            source
        } else {
            golden.join(name)
        };
        fs::copy(&source, destination.join(name)).unwrap();
    }
}

fn assert_text_eq(left: impl AsRef<Path>, right: impl AsRef<Path>) {
    let left = fs::read_to_string(left).unwrap();
    let right = fs::read_to_string(right).unwrap();
    assert_eq!(left, right);
}

#[test]
fn grabpars_returns_expected_value() {
    let tmp = tempfile::tempdir().unwrap();
    copy_files(&["par.dates"], tmp.path());
    let mut command = Command::cargo_bin("grabpars").unwrap();
    command
        .current_dir(tmp.path())
        .args(["-p", "par.dates", "-x", "output:"])
        .assert()
        .success()
        .stdout("Toy.out\n");
}

#[test]
fn dowtjack_matches_golden_output() {
    let tmp = tempfile::tempdir().unwrap();
    copy_files(&["Toy.jin"], tmp.path());
    Command::cargo_bin("dowtjack")
        .unwrap()
        .current_dir(tmp.path())
        .args(["-i", "Toy.jin", "-o", "Toy.jout", "-m", "0.0"])
        .assert()
        .success();
    assert_text_eq(tmp.path().join("Toy.jout"), golden_dir().join("Toy.jout"));
}

#[test]
fn simpjack2_prints_summary() {
    let tmp = tempfile::tempdir().unwrap();
    copy_files(&["Toy.jin"], tmp.path());
    Command::cargo_bin("simpjack2")
        .unwrap()
        .current_dir(tmp.path())
        .args(["-i", "Toy.jin", "-m", "0.0"])
        .assert()
        .success();
}

#[test]
fn dates_expfit_matches_golden_fit() {
    let tmp = tempfile::tempdir().unwrap();
    copy_files(&["Toy.out"], tmp.path());
    Command::cargo_bin("dates_expfit")
        .unwrap()
        .current_dir(tmp.path())
        .args([
            "-i", "Toy.out", "-n", "1", "-o", "Toy.fit", "-l", "0.45", "-s", "0.005", "-a", "-c",
            "3", "-r", "77",
        ])
        .assert()
        .success();
    assert_text_eq(tmp.path().join("Toy.fit"), golden_dir().join("Toy.fit"));
}

#[test]
fn dates_plot_writes_expected_artifacts() {
    let tmp = tempfile::tempdir().unwrap();
    copy_files(&["Toy.out"], tmp.path());
    Command::cargo_bin("dates_plot")
        .unwrap()
        .current_dir(tmp.path())
        .args([
            "-i", "Toy", "-s", "0.005", "-l", "0.45", "-h", "20", "-r", "77", "-a", "-c", "3",
        ])
        .assert()
        .success();
    assert_text_eq(tmp.path().join("Toy.fit"), golden_dir().join("Toy.fit"));
    assert_text_eq(tmp.path().join("Toy.xtxt"), golden_dir().join("Toy.xtxt"));
    assert!(tmp.path().join("Toy.ps").exists());
    assert!(tmp.path().join("Toy.pdf").exists());
}

#[test]
fn run_dates_expfit_matches_golden_output() {
    let tmp = tempfile::tempdir().unwrap();
    copy_files(&["par.dates", "Toy.out"], tmp.path());
    Command::cargo_bin("run_dates_expfit")
        .unwrap()
        .current_dir(tmp.path())
        .args(["-p", "par.dates"])
        .assert()
        .success();
    assert_text_eq(
        tmp.path().join("Toy:expfit.out"),
        golden_dir().join("Toy:expfit.out"),
    );
}

#[test]
fn dates_jackknife_matches_golden_text_outputs() {
    let tmp = tempfile::tempdir().unwrap();
    copy_files(
        &["par.dates", "toy.snp", "Toy.out", "Toy.out:1", "Toy.out:2"],
        tmp.path(),
    );
    Command::cargo_bin("dates_jackknife")
        .unwrap()
        .current_dir(tmp.path())
        .args(["-p", "par.dates", "-m", "toy.snp", "-r", "77", "-a"])
        .assert()
        .success();
    assert_text_eq(tmp.path().join("Toy.jin"), golden_dir().join("Toy.jin"));
    assert_text_eq(tmp.path().join("Toy.jout"), golden_dir().join("Toy.jout"));
}

#[test]
fn dates_matches_golden_outputs() {
    let tmp = tempfile::tempdir().unwrap();
    copy_files(
        &["par.dates", "toy.geno", "toy.snp", "toy.ind", "poplist.txt"],
        tmp.path(),
    );
    Command::cargo_bin("dates")
        .unwrap()
        .current_dir(tmp.path())
        .args(["-p", "par.dates"])
        .assert()
        .success();
    for name in [
        "Toy.out",
        "Toy.out:1",
        "Toy.out:2",
        "Toy.jin",
        "Toy.jout",
        "Toy.fit",
        "Toy.xtxt",
        "Mix:log",
    ] {
        assert_text_eq(tmp.path().join(name), golden_dir().join(name));
    }
    assert!(tmp.path().join("Toy.ps").exists());
    assert!(tmp.path().join("Toy.pdf").exists());
}
