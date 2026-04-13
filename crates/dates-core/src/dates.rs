//! Main DATES runtime.
//!
//! This module ports the supported execution path from `src/dates.c`,
//! retaining the job model, covariance output layout, per-chromosome jackknife
//! files, and helper-workflow integration. The current Rust implementation
//! targets text Eigenstrat inputs, which is sufficient for the Rust fixture
//! corpus and the documented self-contained workflow in this workspace.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow, bail};

use crate::config::DatesParams;
use crate::corr::Corr;
use crate::dataset::{AdmixJob, Dataset, load_individual_values, load_weights};
use crate::fft::{auto_positive, cross_positive};
use crate::workflow::{DatesJackknifeRequest, run_dates_expfit_from_par, run_dates_jackknife};

#[derive(Clone, Debug)]
struct SelectedSnp {
    snp_index: usize,
    chrom: i32,
    genpos: f64,
    qbin_tag: i32,
    parent_a_freq: f64,
    parent_b_freq: f64,
    weight: f64,
}

#[derive(Clone, Copy, Debug)]
struct PresentSnp {
    chrom: i32,
    genpos: f64,
    qbin_tag: usize,
    genotype: f64,
    parent_a_freq: f64,
    parent_b_freq: f64,
    weight: f64,
}

/// Run the Rust `dates` entrypoint.
pub fn run_dates(par_path: &Path, verbose: bool) -> Result<()> {
    println!();
    println!("## DATES.  Version 753");
    println!();

    let par_path = par_path
        .canonicalize()
        .unwrap_or_else(|_| par_path.to_path_buf());
    let params = DatesParams::load(&par_path)?;
    if verbose {
        eprintln!("parameter file: {}", par_path.display());
    }

    let badsnp = params.badsnpname.as_ref().map(Path::new);
    let dataset = Dataset::load(
        resolve_param_path(&par_path, &params.genotypename),
        resolve_param_path(&par_path, &params.snpname),
        resolve_param_path(&par_path, &params.indivname),
        badsnp,
    )?;
    if params.checkmap && !dataset.has_real_map() {
        bail!("running DATES without a real map; set checkmap: NO if that is intentional");
    }
    let jobs = load_jobs_resolved(&params, &par_path)?;
    let weight_map = if let Some(path) = &params.weightname {
        load_weights(resolve_param_path(&par_path, path))?
    } else {
        BTreeMap::new().into_iter().collect()
    };
    let timeoffset_map = if let Some(path) = &params.timeoffsetname {
        load_individual_values(resolve_param_path(&par_path, path))?
    } else {
        BTreeMap::new().into_iter().collect()
    };
    for job in &jobs {
        run_job(
            &par_path,
            &dataset,
            &params,
            job,
            &weight_map,
            &timeoffset_map,
            verbose,
        )?;
    }
    Ok(())
}

fn run_job(
    par_path: &Path,
    dataset: &Dataset,
    params: &DatesParams,
    job: &AdmixJob,
    weight_map: &std::collections::HashMap<String, f64>,
    timeoffset_map: &std::collections::HashMap<String, f64>,
    verbose: bool,
) -> Result<()> {
    let output_dir = if let Some(dir) = &job.output_dir {
        fs::create_dir_all(dir).with_context(|| format!("failed to create {}", dir.display()))?;
        if let Some(file_name) = par_path.file_name() {
            let target = dir.join(file_name);
            let _ = fs::copy(par_path, &target);
        }
        dir.clone()
    } else {
        std::env::current_dir()?
    };

    let prefix = if job.output_dir.is_none() {
        params
            .output
            .as_deref()
            .map(strip_out_suffix)
            .map(ToOwned::to_owned)
            .unwrap_or_else(|| job.admixed.clone())
    } else {
        job.admixed.clone()
    };
    let output_prefix = output_dir.join(&prefix);
    let main_output = output_prefix.with_extension("out");
    let logfit = output_dir.join(format!("{}:log", job.admixed));

    if verbose {
        eprintln!(
            "running job {} {} {} -> {}",
            job.source_a,
            job.source_b,
            job.admixed,
            output_prefix.display()
        );
    }

    let parent_a = dataset
        .individuals
        .iter()
        .enumerate()
        .filter_map(|(index, individual)| (individual.egroup == job.source_a).then_some(index))
        .collect::<Vec<_>>();
    let parent_b = dataset
        .individuals
        .iter()
        .enumerate()
        .filter_map(|(index, individual)| (individual.egroup == job.source_b).then_some(index))
        .collect::<Vec<_>>();
    let admixed = dataset
        .individuals
        .iter()
        .enumerate()
        .filter_map(|(index, individual)| (individual.egroup == job.admixed).then_some(index))
        .collect::<Vec<_>>();

    if parent_a.is_empty() || parent_b.is_empty() {
        bail!("admixing population has no samples");
    }
    if admixed.is_empty() {
        bail!("no admixed samples found");
    }

    let mut selected = Vec::new();
    for (snp_index, snp) in dataset.snps.iter().enumerate() {
        if snp.chrom < 1 || snp.chrom > 22 {
            continue;
        }
        if let Some(chrom) = params.chrom
            && snp.chrom != chrom
        {
            continue;
        }
        if let Some(nochrom) = params.nochrom
            && snp.chrom == nochrom
        {
            continue;
        }
        let parent_a_freq = allele_frequency(snp, &parent_a)?;
        let parent_b_freq = allele_frequency(snp, &parent_b)?;
        let diff = parent_a_freq - parent_b_freq;
        if diff.abs() < 0.001 {
            continue;
        }
        let mean_freq = 0.5 * (parent_a_freq + parent_b_freq);
        let denom = mean_freq * (1.0 - mean_freq);
        if denom < 0.001 {
            continue;
        }
        let mut weight = match params.runmode {
            0 => diff / mean_freq.sqrt(),
            1 => diff,
            2 => 1.0,
            other => bail!("unsupported runmode {other}"),
        };
        if let Some(override_weight) = weight_map.get(&snp.id) {
            weight = *override_weight;
        }
        selected.push(SelectedSnp {
            snp_index,
            chrom: snp.chrom,
            genpos: snp.genpos,
            qbin_tag: -1,
            parent_a_freq,
            parent_b_freq,
            weight,
        });
    }
    if selected.is_empty() {
        bail!("zero usable SNPs after filtering");
    }
    if params.qbin > 0 {
        assign_qbins(&mut selected, params.binsize, params.qbin);
    }

    let mut centered_offsets = BTreeMap::new();
    if !timeoffset_map.is_empty() {
        let mean = admixed
            .iter()
            .map(|index| {
                timeoffset_map
                    .get(&dataset.individuals[*index].id)
                    .copied()
                    .unwrap_or(0.0)
            })
            .sum::<f64>()
            / admixed.len() as f64;
        for index in &admixed {
            let individual = &dataset.individuals[*index];
            centered_offsets.insert(
                *index,
                timeoffset_map.get(&individual.id).copied().unwrap_or(0.0) - mean,
            );
        }
    }

    let num_bins = (params.maxdis / params.binsize).round() as usize + 5;
    let mut chrom_corr = vec![vec![Corr::default(); num_bins]; 23];
    if params.qbin > 0 {
        run_qbin_mode(
            dataset,
            &selected,
            &admixed,
            &centered_offsets,
            params,
            &mut chrom_corr,
        )?;
    } else {
        run_direct_mode(
            dataset,
            &selected,
            &admixed,
            &centered_offsets,
            params,
            &mut chrom_corr,
        )?;
    }

    let global = sum_chrom_corr(&chrom_corr, num_bins);
    dump_output(
        &main_output,
        &global,
        num_bins.saturating_sub(5),
        params.binsize,
        &format!(
            " ##Z-score and correlation:: {}  binsize: {:12.6}",
            job.admixed, params.binsize
        ),
        params.runmode,
    )?;
    if params.jackknife {
        for (chrom, chrom_bins) in chrom_corr.iter().enumerate().take(23).skip(1) {
            let mut leave_out = Vec::with_capacity(num_bins);
            for bin in 0..num_bins {
                leave_out.push(global[bin].minus(chrom_bins[bin])?);
            }
            dump_output(
                Path::new(&format!("{}:{}", main_output.display(), chrom)),
                &leave_out,
                num_bins.saturating_sub(5),
                params.binsize,
                &format!("## Jackknife output: chrom {}", chrom),
                params.runmode,
            )?;
        }
    }

    if params.runfit {
        let cwd = std::env::current_dir()?;
        std::env::set_current_dir(&output_dir)?;
        let summary = if params.jackknife {
            run_dates_jackknife(&DatesJackknifeRequest {
                par_path: par_path.to_path_buf(),
                data_col: 3,
                low_cm: params.lovalfit,
                high_cm: 20.0,
                snp_override: Some(resolve_param_path(par_path, &params.snpname)),
                admix_override: Some(job.admixed.clone()),
                affine: params.afffit,
                seed: params.seed,
            })?
        } else {
            let _ = run_dates_expfit_from_par(
                par_path,
                3,
                params.lovalfit,
                params.afffit,
                params.seed,
                Some(&job.admixed),
            )?;
            crate::dataset::JackknifeSummary {
                mean: 0.0,
                std_err: 0.0,
            }
        };
        std::env::set_current_dir(cwd)?;
        let mut log_body = String::new();
        log_body.push_str("calling fit!...\n");
        log_body.push_str(&format!("prefix: {}\n", prefix));
        if params.jackknife {
            log_body.push_str(&format!(
                "jackknife summary: {:9.3} {:9.3}\n",
                summary.mean, summary.std_err
            ));
        }
        fs::write(&logfit, log_body)
            .with_context(|| format!("failed to write {}", logfit.display()))?;
    }
    Ok(())
}

fn run_direct_mode(
    dataset: &Dataset,
    selected: &[SelectedSnp],
    admixed: &[usize],
    centered_offsets: &BTreeMap<usize, f64>,
    params: &DatesParams,
    chrom_corr: &mut [Vec<Corr>],
) -> Result<()> {
    for individual_index in admixed {
        let timeoffset = centered_offsets
            .get(individual_index)
            .copied()
            .unwrap_or(0.0);
        let present = build_present_snps(dataset, selected, *individual_index)?;
        if present.is_empty() {
            continue;
        }
        let values = build_weighted_values(&present)?;
        for left in 0..present.len() {
            for right in left + 1..present.len() {
                if present[left].chrom != present[right].chrom {
                    break;
                }
                let dis = present[right].genpos - present[left].genpos;
                if dis >= params.maxdis {
                    break;
                }
                let bin = (dis / params.binsize) as usize;
                let y1 = values[left] * (-timeoffset * dis).exp();
                let y2 = values[right];
                chrom_corr[present[left].chrom as usize][bin].add(y1, y2);
            }
        }
    }
    Ok(())
}

fn run_qbin_mode(
    dataset: &Dataset,
    selected: &[SelectedSnp],
    admixed: &[usize],
    centered_offsets: &BTreeMap<usize, f64>,
    params: &DatesParams,
    chrom_corr: &mut [Vec<Corr>],
) -> Result<()> {
    let num_bins = chrom_corr[1].len();
    let num_qbins = selected
        .last()
        .map(|snp| snp.qbin_tag.max(0) as usize + 1)
        .unwrap_or(0);
    let num_dbins = num_bins * params.qbin;
    let diff_max = ((params.qbin as f64) * params.maxdis / params.binsize).round() as usize;
    let mut ddcbins = vec![vec![vec![0.0; num_dbins]; 7]; 23];
    for individual_index in admixed {
        let timeoffset = centered_offsets
            .get(individual_index)
            .copied()
            .unwrap_or(0.0);
        let present = build_present_snps(dataset, selected, *individual_index)?;
        if present.is_empty() {
            continue;
        }
        let values = build_weighted_values(&present)?;
        let mut z0 = vec![0.0; num_qbins];
        let mut z1 = vec![0.0; num_qbins];
        let mut z2 = vec![0.0; num_qbins];
        let mut ranges = BTreeMap::<i32, (usize, usize)>::new();
        for (index, snp) in present.iter().enumerate() {
            let tag = snp.qbin_tag;
            let weighted = values[index];
            z0[tag] += 1.0;
            z1[tag] += weighted;
            z2[tag] += weighted * weighted;
            let entry = ranges.entry(snp.chrom).or_insert((tag, tag));
            entry.0 = entry.0.min(tag);
            entry.1 = entry.1.max(tag);
        }
        for (chrom, (start, end)) in ranges {
            if end <= start {
                continue;
            }
            let z0s = &z0[start..=end];
            let z1s = &z1[start..=end];
            let z2s = &z2[start..=end];
            let dd00 = pad_corr(auto_positive(z0s, diff_max), diff_max + 1);
            let dd01 = pad_corr(cross_positive(z0s, z1s, diff_max), diff_max + 1);
            let dd10 = pad_corr(cross_positive(z1s, z0s, diff_max), diff_max + 1);
            let dd11 = pad_corr(auto_positive(z1s, diff_max), diff_max + 1);
            let dd02 = pad_corr(cross_positive(z0s, z2s, diff_max), diff_max + 1);
            let dd20 = pad_corr(cross_positive(z2s, z0s, diff_max), diff_max + 1);
            for d in 0..=diff_max.min(num_dbins.saturating_sub(1)) {
                ddcbins[chrom as usize][0][d] += dd00[d];
                ddcbins[chrom as usize][1][d] += dd01[d];
                ddcbins[chrom as usize][2][d] += dd02[d];
                ddcbins[chrom as usize][3][d] += dd10[d];
                ddcbins[chrom as usize][4][d] += dd11[d];
                ddcbins[chrom as usize][6][d] += dd20[d];
            }
        }
        let _ = timeoffset;
    }
    let dbinsize = params.binsize / params.qbin as f64;
    for (chrom, chrom_bins) in ddcbins.iter().enumerate().take(23).skip(1) {
        for (d, dd00) in chrom_bins[0].iter().enumerate().take(num_dbins).skip(1) {
            let dd00 = *dd00;
            if dd00 < 0.5 {
                continue;
            }
            let ys = d as f64 * dbinsize;
            let bin = (ys / params.binsize) as usize;
            if bin >= chrom_corr[chrom].len() {
                continue;
            }
            chrom_corr[chrom][bin].add_summary(
                dd00,
                chrom_bins[1][d],
                chrom_bins[3][d],
                chrom_bins[4][d],
                chrom_bins[2][d],
                chrom_bins[6][d],
            );
        }
    }
    Ok(())
}

fn build_present_snps(
    dataset: &Dataset,
    selected: &[SelectedSnp],
    individual_index: usize,
) -> Result<Vec<PresentSnp>> {
    let mut present = Vec::new();
    for snp in selected {
        let genotype = dataset.snps[snp.snp_index]
            .gtypes
            .get(individual_index)
            .copied()
            .ok_or_else(|| anyhow!("missing genotype"))?;
        if genotype < 0 {
            continue;
        }
        present.push(PresentSnp {
            chrom: snp.chrom,
            genpos: snp.genpos,
            qbin_tag: snp.qbin_tag.max(0) as usize,
            genotype: genotype as f64 / 2.0,
            parent_a_freq: snp.parent_a_freq,
            parent_b_freq: snp.parent_b_freq,
            weight: snp.weight,
        });
    }
    Ok(present)
}

fn build_weighted_values(present: &[PresentSnp]) -> Result<Vec<f64>> {
    let w0 = present.iter().map(|snp| snp.genotype).collect::<Vec<_>>();
    let w1 = present
        .iter()
        .map(|snp| snp.parent_a_freq)
        .collect::<Vec<_>>();
    let w2 = present
        .iter()
        .map(|snp| snp.parent_b_freq)
        .collect::<Vec<_>>();
    let ww1 = w0
        .iter()
        .zip(w2.iter())
        .map(|(lhs, rhs)| lhs - rhs)
        .collect::<Vec<_>>();
    let ww2 = w1
        .iter()
        .zip(w2.iter())
        .map(|(lhs, rhs)| lhs - rhs)
        .collect::<Vec<_>>();
    let denom = dot(&ww2, &ww2);
    if denom <= 1.0e-20 {
        bail!("invalid ancestry-fit denominator");
    }
    let coeff = dot(&ww1, &ww2) / denom;
    Ok(present
        .iter()
        .zip(w1.iter().zip(w2.iter()))
        .map(|(snp, (freq_a, freq_b))| {
            let prediction = coeff * freq_a + (1.0 - coeff) * freq_b;
            (snp.genotype - prediction) * snp.weight
        })
        .collect())
}

fn allele_frequency(snp: &crate::dataset::Snp, individuals: &[usize]) -> Result<f64> {
    let mut counts = [0.0; 3];
    for index in individuals {
        let genotype = *snp
            .gtypes
            .get(*index)
            .ok_or_else(|| anyhow!("missing genotype"))?;
        if genotype < 0 {
            continue;
        }
        counts[genotype as usize] += 1.0;
    }
    let total = counts.iter().sum::<f64>() * 2.0;
    if total <= 0.01 {
        bail!("no parental genotypes");
    }
    Ok((2.0 * counts[2] + counts[1]) / total)
}

fn assign_qbins(selected: &mut [SelectedSnp], binsize: f64, qbin: usize) {
    if selected.is_empty() {
        return;
    }
    let qb = binsize / qbin as f64;
    let mut ydis = 0.0;
    let mut lastgenpos = selected[0].genpos;
    let mut chrom = selected[0].chrom;
    selected[0].qbin_tag = 0;
    for snp in selected.iter_mut().skip(1) {
        if snp.chrom != chrom {
            ydis += 5.0;
            chrom = snp.chrom;
        } else {
            ydis += snp.genpos - lastgenpos;
        }
        lastgenpos = snp.genpos;
        snp.qbin_tag = (ydis / qb).floor() as i32;
    }
}

fn sum_chrom_corr(chrom_corr: &[Vec<Corr>], num_bins: usize) -> Vec<Corr> {
    let mut global = vec![Corr::default(); num_bins];
    for chrom_bins in chrom_corr.iter().take(23).skip(1) {
        for bin in 0..num_bins {
            global[bin] = global[bin].plus(chrom_bins[bin]);
        }
    }
    global
}

fn dump_output(
    path: &Path,
    bins: &[Corr],
    len: usize,
    binsize: f64,
    header: &str,
    runmode: i32,
) -> Result<()> {
    let mut out = String::new();
    out.push_str(header);
    out.push('\n');
    for (index, corr) in bins.iter().take(len).enumerate() {
        let mut corr = *corr;
        out.push_str(&format!("{:9.3} ", 100.0 * (index + 1) as f64 * binsize));
        if runmode != 2 {
            let _ = corr.calculate(true, false)?;
            let regression = corr.s12 / (corr.s11 + 1.0e-20);
            out.push_str(&format!(
                "{:15.9} {:12.6} {:12.6} {:12.0} ",
                corr.v12, regression, corr.corr, corr.s0
            ));
        }
        out.push('\n');
    }
    fs::write(path, out).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right.iter())
        .map(|(lhs, rhs)| lhs * rhs)
        .sum()
}

fn pad_corr(mut values: Vec<f64>, target_len: usize) -> Vec<f64> {
    if values.len() < target_len {
        values.resize(target_len, 0.0);
    }
    values
}

fn strip_out_suffix(value: &str) -> &str {
    value.strip_suffix(".out").unwrap_or(value)
}

fn resolve_param_path(par_path: &Path, raw: &str) -> PathBuf {
    let path = Path::new(raw);
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        par_path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(path)
    }
}

fn load_jobs_resolved(params: &DatesParams, par_path: &Path) -> Result<Vec<AdmixJob>> {
    if let Some(admixlist) = &params.admixlist {
        let path = resolve_param_path(par_path, admixlist);
        let raw = fs::read_to_string(&path)
            .with_context(|| format!("failed to read admix list {}", path.display()))?;
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
                output_dir: Some(resolve_param_path(par_path, parts[3])),
            });
        }
        return Ok(jobs);
    }
    let poplist = params
        .poplistname
        .as_deref()
        .ok_or_else(|| anyhow!("missing poplistname"))?;
    let path = resolve_param_path(par_path, poplist);
    let raw =
        fs::read_to_string(&path).with_context(|| format!("failed to read {}", path.display()))?;
    let groups = raw
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(ToOwned::to_owned)
        .collect::<Vec<_>>();
    if groups.len() != 2 {
        bail!("poplistname must contain exactly two groups");
    }
    Ok(vec![AdmixJob {
        source_a: groups[0].clone(),
        source_b: groups[1].clone(),
        admixed: params
            .admixpop
            .clone()
            .ok_or_else(|| anyhow!("missing admixpop"))?,
        output_dir: None,
    }])
}
