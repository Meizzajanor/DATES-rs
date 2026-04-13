//! Dataset loading and typed data model for DATES.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow, bail};

use crate::config::DatesParams;

/// One individual entry from an Eigenstrat `.ind` file.
#[derive(Clone, Debug)]
pub struct Individual {
    pub id: String,
    pub gender: char,
    pub egroup: String,
    pub qval: f64,
}

/// One SNP entry and its loaded genotypes.
#[derive(Clone, Debug)]
pub struct Snp {
    pub id: String,
    pub chrom: i32,
    pub genpos: f64,
    pub physpos: f64,
    pub alleles: [char; 2],
    pub gtypes: Vec<i8>,
    pub qbin_tag: i32,
}

/// A loaded DATES dataset.
#[derive(Clone, Debug)]
pub struct Dataset {
    pub individuals: Vec<Individual>,
    pub snps: Vec<Snp>,
}

/// One admix job entry from `admixlist`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AdmixJob {
    pub source_a: String,
    pub source_b: String,
    pub admixed: String,
    pub output_dir: Option<PathBuf>,
}

/// One written covariance row.
#[derive(Clone, Debug, PartialEq)]
pub struct CovarianceBin {
    pub distance_cm: f64,
    pub covariance: f64,
    pub regression: f64,
    pub correlation: f64,
    pub pairs: f64,
}

/// Fitting request used by helper workflows.
#[derive(Clone, Debug)]
pub struct FitRequest {
    pub input: PathBuf,
    pub output: Option<PathBuf>,
    pub num_exp: usize,
    pub data_col: usize,
    pub low_cm: f64,
    pub high_cm: f64,
    pub step_morgans: Option<f64>,
    pub add_x: f64,
    pub affine: bool,
    pub seed: u64,
}

/// One fitted row emitted by `dates_expfit`.
#[derive(Clone, Debug, PartialEq)]
pub struct FitRow {
    pub distance_cm: f64,
    pub observed: f64,
    pub fitted: f64,
    pub residual: f64,
}

/// Result returned by the exponential fitter.
#[derive(Clone, Debug)]
pub struct FitResult {
    pub rows: Vec<FitRow>,
    pub error_sd: f64,
    pub halflives: Vec<f64>,
    pub mean_generations: Vec<f64>,
    pub basis: Vec<f64>,
    pub coefficients: Vec<f64>,
    pub step_morgans: f64,
}

/// One jackknife sample row.
#[derive(Clone, Debug, PartialEq)]
pub struct JackknifeSample {
    pub block: i32,
    pub weight: f64,
    pub estimate: f64,
}

/// Weighted jackknife summary.
#[derive(Clone, Debug, PartialEq)]
pub struct JackknifeSummary {
    pub mean: f64,
    pub std_err: f64,
}

impl Dataset {
    /// Load a text-Eigenstrat dataset for the DATES runtime.
    pub fn load(
        genotype_path: impl AsRef<Path>,
        snp_path: impl AsRef<Path>,
        indiv_path: impl AsRef<Path>,
        badsnp_path: Option<&Path>,
    ) -> Result<Self> {
        let individuals = load_individuals(indiv_path.as_ref())?;
        let mut snps = load_snps(snp_path.as_ref(), badsnp_path)?;
        load_text_genotypes(genotype_path.as_ref(), &mut snps, individuals.len())?;
        Ok(Self { individuals, snps })
    }

    /// Return true when the genetic map does not look like a physical-position placeholder.
    pub fn has_real_map(&self) -> bool {
        self.snps.iter().take(10).any(|snp| {
            let physical = snp.physpos / 1.0e8;
            (snp.genpos - physical).abs() > 0.001
        })
    }
}

/// Load `admixlist` or `poplistname` from a parsed `dates` parameter set.
pub fn load_jobs(params: &DatesParams) -> Result<Vec<AdmixJob>> {
    if let Some(admixlist) = &params.admixlist {
        let raw = fs::read_to_string(admixlist)
            .with_context(|| format!("failed to read admix list {}", admixlist))?;
        let mut jobs = Vec::new();
        for line in raw.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let parts: Vec<_> = trimmed.split_whitespace().collect();
            if parts.len() != 4 {
                bail!("admixlist line must have 4 fields: {trimmed}");
            }
            jobs.push(AdmixJob {
                source_a: parts[0].to_owned(),
                source_b: parts[1].to_owned(),
                admixed: parts[2].to_owned(),
                output_dir: Some(PathBuf::from(parts[3])),
            });
        }
        return Ok(jobs);
    }
    let admixed = params
        .admixpop
        .clone()
        .ok_or_else(|| anyhow!("missing admixpop: for single-job execution"))?;
    let poplist = params
        .poplistname
        .clone()
        .ok_or_else(|| anyhow!("missing poplistname: for single-job execution"))?;
    let raw = fs::read_to_string(&poplist).with_context(|| format!("failed to read {poplist}"))?;
    let groups: Vec<_> = raw
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(ToOwned::to_owned)
        .collect();
    if groups.len() != 2 {
        bail!("poplistname must contain exactly two groups");
    }
    Ok(vec![AdmixJob {
        source_a: groups[0].clone(),
        source_b: groups[1].clone(),
        admixed,
        output_dir: None,
    }])
}

/// Count SNPs per chromosome from an `.snp` file.
pub fn chromosome_counts(path: impl AsRef<Path>) -> Result<BTreeMap<i32, usize>> {
    let mut counts = BTreeMap::new();
    for line in fs::read_to_string(path.as_ref())?.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<_> = trimmed.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        let chrom = parts[1]
            .parse::<i32>()
            .with_context(|| format!("invalid chromosome in {trimmed}"))?;
        *counts.entry(chrom).or_insert(0) += 1;
    }
    Ok(counts)
}

/// Load optional per-individual numeric values by individual id.
pub fn load_individual_values(path: impl AsRef<Path>) -> Result<HashMap<String, f64>> {
    let mut values = HashMap::new();
    for line in fs::read_to_string(path.as_ref())?.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<_> = trimmed.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        values.insert(parts[0].to_owned(), parts[1].parse::<f64>()?);
    }
    Ok(values)
}

/// Load optional per-SNP weights by SNP id.
pub fn load_weights(path: impl AsRef<Path>) -> Result<HashMap<String, f64>> {
    let mut values = HashMap::new();
    for line in fs::read_to_string(path.as_ref())?.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<_> = trimmed.split_whitespace().collect();
        if parts.len() < 2 {
            continue;
        }
        values.insert(parts[0].to_owned(), parts[1].parse::<f64>()?);
    }
    Ok(values)
}

fn load_individuals(path: &Path) -> Result<Vec<Individual>> {
    let raw =
        fs::read_to_string(path).with_context(|| format!("failed to read {}", path.display()))?;
    let mut individuals = Vec::new();
    for line in raw.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<_> = trimmed.split_whitespace().collect();
        if parts.len() < 3 {
            bail!("invalid indiv row: {trimmed}");
        }
        individuals.push(Individual {
            id: parts[0].to_owned(),
            gender: parts[1].chars().next().unwrap_or('U'),
            egroup: parts[2].to_owned(),
            qval: 0.0,
        });
    }
    Ok(individuals)
}

fn load_snps(path: &Path, badsnp_path: Option<&Path>) -> Result<Vec<Snp>> {
    let bad_snps: BTreeSet<String> = if let Some(path) = badsnp_path {
        fs::read_to_string(path)?
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .map(ToOwned::to_owned)
            .collect()
    } else {
        BTreeSet::new()
    };
    let raw =
        fs::read_to_string(path).with_context(|| format!("failed to read {}", path.display()))?;
    let mut snps = Vec::new();
    for line in raw.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<_> = trimmed.split_whitespace().collect();
        if parts.len() < 4 {
            bail!("invalid snp row: {trimmed}");
        }
        if bad_snps.contains(parts[0]) {
            continue;
        }
        let alleles = if parts.len() >= 6 {
            [
                parts[4].chars().next().unwrap_or('X'),
                parts[5].chars().next().unwrap_or('X'),
            ]
        } else {
            ['X', 'X']
        };
        snps.push(Snp {
            id: parts[0].to_owned(),
            chrom: parts[1].parse::<i32>()?,
            genpos: parts[2].parse::<f64>()?,
            physpos: parts[3].parse::<f64>()?,
            alleles,
            gtypes: Vec::new(),
            qbin_tag: -1,
        });
    }
    Ok(snps)
}

fn load_text_genotypes(path: &Path, snps: &mut [Snp], num_individuals: usize) -> Result<()> {
    let file =
        File::open(path).with_context(|| format!("failed to read {}", path.display()))?;
    let mut reader = BufReader::new(file);

    // Peek at the first four bytes to detect packed Eigenstrat format.
    let buf = reader.fill_buf()?;
    if buf.starts_with(b"GENO") {
        bail!(
            "packed Eigenstrat input is not yet supported in DATES-rs: {}",
            path.display()
        );
    }

    let mut line_buf = String::new();
    let mut line_count: usize = 0;
    for snp in snps.iter_mut() {
        line_buf.clear();
        let bytes_read = reader
            .read_line(&mut line_buf)
            .with_context(|| format!("failed to read line from geno file {}", path.display()))?;
        if bytes_read == 0 {
            bail!(
                "geno line count {} does not match snp count {} for {}",
                line_count,
                snps.len(),
                path.display()
            );
        }
        line_count += 1;
        let trimmed = line_buf.trim();
        if trimmed.len() != num_individuals {
            bail!(
                "geno row length {} does not match indiv count {} for SNP {}",
                trimmed.len(),
                num_individuals,
                snp.id
            );
        }
        let mut gtypes = Vec::with_capacity(num_individuals);
        for byte in trimmed.bytes() {
            let value = match byte {
                b'0' => 0,
                b'1' => 1,
                b'2' => 2,
                b'9' => -1,
                other => bail!("invalid genotype byte {}", other as char),
            };
            gtypes.push(value);
        }
        snp.gtypes = gtypes;
    }
    // Verify that the file has no extra lines beyond what was expected.
    line_buf.clear();
    if reader.read_line(&mut line_buf)? > 0 {
        bail!(
            "geno line count {} does not match snp count {} for {}",
            line_count + 1,
            snps.len(),
            path.display()
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_admix_jobs() {
        let params = DatesParams {
            par_path: PathBuf::from("x"),
            genotypename: "g".into(),
            snpname: "s".into(),
            indivname: "i".into(),
            poplistname: Some("unused".into()),
            weightname: None,
            timeoffsetname: None,
            admixpop: Some("SIM".into()),
            admixlist: None,
            badsnpname: None,
            output: None,
            runmode: 0,
            ldmode: false,
            seed: 0,
            nxlim: 0,
            minparentcount: 0,
            norun: false,
            flatmode: false,
            zdipcorrmode: false,
            jackknife: true,
            maxdis: 0.1,
            binsize: 0.01,
            chithresh: 0.0,
            checkmap: true,
            qbin: 0,
            runfit: true,
            afffit: true,
            lovalfit: 0.45,
            chrom: None,
            nochrom: None,
            numchrom: 22,
        };
        assert_eq!(params.admixpop.as_deref(), Some("SIM"));
    }
}
