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

fn write_chrom23_fixture(destination: &Path) {
    fs::write(destination.join("toy.ind"), "A U SRC1\nB U SRC2\nM U MIX\n").unwrap();
    fs::write(destination.join("poplist.txt"), "SRC1\nSRC2\n").unwrap();
    fs::write(
        destination.join("toy.snp"),
        "rs1 23 0.000 100 A G\nrs2 23 0.005 200 A G\n",
    )
    .unwrap();
    fs::write(destination.join("toy.geno"), "201\n201\n").unwrap();
    fs::write(
        destination.join("par.dates"),
        "genotypename: toy.geno\n\
         snpname: toy.snp\n\
         indivname: toy.ind\n\
         poplistname: poplist.txt\n\
         admixpop: MIX\n\
         output: Chrom23.out\n\
         binsize: 0.005\n\
         maxdis: 0.010\n\
         seed: 77\n\
         runmode: 1\n\
         checkmap: NO\n\
         numchrom: 23\n\
         qbin: 0\n\
         jackknife: NO\n\
         runfit: NO\n",
    )
    .unwrap();
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
fn run_dates_expfit_reads_par_from_another_directory() {
    let par_dir = tempfile::tempdir().unwrap();
    let work_dir = tempfile::tempdir().unwrap();
    copy_files(&["par.dates"], par_dir.path());
    copy_files(&["Toy.out"], work_dir.path());
    Command::cargo_bin("run_dates_expfit")
        .unwrap()
        .current_dir(work_dir.path())
        .args(["-p", par_dir.path().join("par.dates").to_str().unwrap()])
        .assert()
        .success();
    assert_text_eq(
        work_dir.path().join("Toy:expfit.out"),
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
fn dates_jackknife_resolves_parameter_paths_relative_to_par_file() {
    let par_dir = tempfile::tempdir().unwrap();
    let work_dir = tempfile::tempdir().unwrap();
    copy_files(&["par.dates", "toy.snp"], par_dir.path());
    copy_files(&["Toy.out", "Toy.out:1", "Toy.out:2"], work_dir.path());
    Command::cargo_bin("dates_jackknife")
        .unwrap()
        .current_dir(work_dir.path())
        .args([
            "-p",
            par_dir.path().join("par.dates").to_str().unwrap(),
            "-r",
            "77",
            "-a",
        ])
        .assert()
        .success();
    assert_text_eq(
        work_dir.path().join("Toy.jin"),
        golden_dir().join("Toy.jin"),
    );
    assert_text_eq(
        work_dir.path().join("Toy.jout"),
        golden_dir().join("Toy.jout"),
    );
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

#[test]
fn dates_plot_honors_non_default_x_range() {
    let tmp = tempfile::tempdir().unwrap();
    copy_files(&["Toy.out"], tmp.path());
    Command::cargo_bin("dates_plot")
        .unwrap()
        .current_dir(tmp.path())
        .args([
            "-i", "Toy", "-s", "0.005", "-l", "1", "-h", "5", "-r", "77", "-a", "-c", "3",
        ])
        .assert()
        .success();
    let xtxt = fs::read_to_string(tmp.path().join("Toy.xtxt")).unwrap();
    assert!(xtxt.contains("set xrange [1:5]"));
}

#[test]
fn dates_aggregates_chromosomes_beyond_22() {
    let tmp = tempfile::tempdir().unwrap();
    write_chrom23_fixture(tmp.path());
    Command::cargo_bin("dates")
        .unwrap()
        .current_dir(tmp.path())
        .args(["-p", "par.dates"])
        .assert()
        .success();
    let output = fs::read_to_string(tmp.path().join("Chrom23.out")).unwrap();
    let rows = output.lines().skip(1).collect::<Vec<_>>();
    assert_eq!(rows.len(), 2);
    let first = rows[0].split_whitespace().collect::<Vec<_>>();
    let second = rows[1].split_whitespace().collect::<Vec<_>>();
    assert_eq!(first.last(), Some(&"0"));
    assert_eq!(second.last(), Some(&"1"));
}

#[test]
fn dates_rejects_runmode_two_explicitly() {
    let tmp = tempfile::tempdir().unwrap();
    copy_files(
        &["toy.geno", "toy.snp", "toy.ind", "poplist.txt"],
        tmp.path(),
    );
    fs::write(
        tmp.path().join("par.dates"),
        "genotypename: toy.geno\n\
         snpname: toy.snp\n\
         indivname: toy.ind\n\
         poplistname: poplist.txt\n\
         admixpop: Mix\n\
         output: Toy.out\n\
         binsize: 0.005\n\
         maxdis: 0.025\n\
         seed: 77\n\
         runmode: 2\n\
         checkmap: NO\n\
         numchrom: 2\n\
         qbin: 0\n\
         jackknife: YES\n\
         runfit: YES\n\
         afffit: YES\n\
         lovalfit: 0.45\n",
    )
    .unwrap();

    let output = Command::cargo_bin("dates")
        .unwrap()
        .current_dir(tmp.path())
        .args(["-p", "par.dates"])
        .output()
        .unwrap();
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("runmode 2 is not yet supported end-to-end in DATES-rs"));
}
