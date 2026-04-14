//! Main DATES runtime.
//!
//! This module ports the supported execution path from `src/dates.c`,
//! retaining the job model, covariance output layout, per-chromosome jackknife
//! files, and helper-workflow integration. The current Rust implementation
//! targets text Eigenstrat inputs, which is sufficient for the Rust fixture
//! corpus and the documented self-contained workflow in this workspace.

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};

use crate::config::{DatesParams, OutputPrefix, resolve_optional_param_path, resolve_param_path};
use crate::corr::Corr;
use crate::dataset::{AdmixJob, Dataset, load_individual_values, load_jobs, load_weights};
use crate::fft::{auto_positive, cross_positive};
use crate::workflow::{
    DatesExpfitRequest, DatesJackknifeRequest, run_dates_expfit_from_par_with_paths,
    run_dates_jackknife,
};

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
    if params.runmode == 2 {
        bail!("runmode 2 is not yet supported end-to-end in DATES-rs");
    }

    let badsnp = resolve_optional_param_path(&par_path, params.badsnpname.as_deref());
    let dataset = Dataset::load(
        resolve_param_path(&par_path, &params.genotypename),
        resolve_param_path(&par_path, &params.snpname),
        resolve_param_path(&par_path, &params.indivname),
        badsnp.as_deref(),
    )?;
    if params.checkmap && !dataset.has_real_map() {
        bail!("running DATES without a real map; set checkmap: NO if that is intentional");
    }
    let jobs = load_jobs(&params)?;
    let weight_map = if let Some(path) = &params.weightname {
        load_weights(resolve_param_path(&par_path, path))?
    } else {
        BTreeMap::new().into_iter().collect()
    };
    let timeoffset_map = if let Some(path) = &params.timeoffsetname {
        if params.qbin > 0 {
            bail!("timeoffsetname is not supported with qbin > 0");
        }
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
    let output_prefix = OutputPrefix::resolve(prefix.clone(), &output_dir);
    let main_output = output_prefix.out_path();
    let logfit = output_dir.join(format!("{}:log", job.admixed));

    if verbose {
        eprintln!(
            "running job {} {} {} -> {}",
            job.source_a,
            job.source_b,
            job.admixed,
            output_prefix.path().display()
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

    if parent_a.is_empty() {
        bail!("source population {} has no samples", job.source_a);
    }
    if parent_b.is_empty() {
        bail!("source population {} has no samples", job.source_b);
    }
    if admixed.is_empty() {
        bail!("no admixed samples found");
    }

    let max_chrom = params.numchrom;
    let num_chroms = max_chrom as usize + 1;
    let mut selected = Vec::new();
    for (snp_index, snp) in dataset.snps.iter().enumerate() {
        if snp.chrom < 1 || snp.chrom > max_chrom {
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
    let mut chrom_corr = vec![vec![Corr::default(); num_bins]; num_chroms];
    if params.qbin > 0 {
        run_qbin_mode(dataset, &selected, &admixed, params, &mut chrom_corr)?;
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
        for (chrom, chrom_bins) in chrom_corr.iter().enumerate().take(num_chroms).skip(1) {
            let mut leave_out = Vec::with_capacity(num_bins);
            for bin in 0..num_bins {
                leave_out.push(global[bin].minus(chrom_bins[bin])?);
            }
            dump_output(
                &output_prefix.chrom_out_path(chrom as i32),
                &leave_out,
                num_bins.saturating_sub(5),
                params.binsize,
                &format!("## Jackknife output: chrom {}", chrom),
                params.runmode,
            )?;
        }
    }

    if params.runfit {
        let summary = if params.jackknife {
            run_dates_jackknife(&DatesJackknifeRequest {
                par_path: par_path.to_path_buf(),
                data_col: 3,
                low_cm: params.lovalfit,
                high_cm: 20.0,
                snp_override: Some(resolve_param_path(par_path, &params.snpname)),
                admix_override: Some(job.admixed.clone()),
                output_dir: Some(output_dir.clone()),
                prefix_override: Some(prefix.clone()),
                affine: params.afffit,
                seed: params.seed,
            })?
        } else {
            let _ = run_dates_expfit_from_par_with_paths(&DatesExpfitRequest {
                par_path,
                data_col: 3,
                low_cm: params.lovalfit,
                affine: params.afffit,
                seed: params.seed,
                admix_override: Some(&job.admixed),
                output_dir: Some(&output_dir),
                prefix_override: Some(&prefix),
            })?;
            crate::dataset::JackknifeSummary {
                mean: 0.0,
                std_err: 0.0,
            }
        };
        let mut log_body = String::new();
        log_body.push_str("calling fit!...\n");
        log_body.push_str(&format!("prefix: {}\n", output_prefix.raw()));
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
    let mut present = Vec::with_capacity(selected.len());
    let mut values = Vec::with_capacity(selected.len());
    for individual_index in admixed {
        let timeoffset = centered_offsets
            .get(individual_index)
            .copied()
            .unwrap_or(0.0);
        build_present_snps(dataset, selected, *individual_index, &mut present)?;
        if present.is_empty() {
            continue;
        }
        build_weighted_values(&present, &mut values)?;
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
    let mut ddcbins = vec![vec![vec![0.0; num_dbins]; 7]; chrom_corr.len()];
    let mut present = Vec::with_capacity(selected.len());
    let mut values = Vec::with_capacity(selected.len());
    for individual_index in admixed {
        build_present_snps(dataset, selected, *individual_index, &mut present)?;
        if present.is_empty() {
            continue;
        }
        build_weighted_values(&present, &mut values)?;
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
    }
    let dbinsize = params.binsize / params.qbin as f64;
    for (chrom, chrom_bins) in ddcbins.iter().enumerate().take(chrom_corr.len()).skip(1) {
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
    present: &mut Vec<PresentSnp>,
) -> Result<()> {
    present.clear();
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
    Ok(())
}

fn build_weighted_values(present: &[PresentSnp], values: &mut Vec<f64>) -> Result<()> {
    values.clear();
    let mut ww1_dot_ww2 = 0.0;
    let mut ww2_dot_ww2 = 0.0;
    for snp in present {
        let diff1 = snp.genotype - snp.parent_b_freq;
        let diff2 = snp.parent_a_freq - snp.parent_b_freq;
        ww1_dot_ww2 += diff1 * diff2;
        ww2_dot_ww2 += diff2 * diff2;
    }
    if ww2_dot_ww2 <= 1.0e-20 {
        bail!("invalid ancestry-fit denominator");
    }
    let coeff = ww1_dot_ww2 / ww2_dot_ww2;
    for snp in present {
        let prediction = coeff * snp.parent_a_freq + (1.0 - coeff) * snp.parent_b_freq;
        values.push((snp.genotype - prediction) * snp.weight);
    }
    Ok(())
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
    for chrom_bins in chrom_corr.iter().skip(1) {
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

fn pad_corr(mut values: Vec<f64>, target_len: usize) -> Vec<f64> {
    if values.len() < target_len {
        values.resize(target_len, 0.0);
    }
    values
}

fn strip_out_suffix(value: &str) -> &str {
    value.strip_suffix(".out").unwrap_or(value)
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use std::sync::{Mutex, OnceLock};

    use super::*;

    fn cwd_lock() -> &'static Mutex<()> {
        static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
        LOCK.get_or_init(|| Mutex::new(()))
    }

    #[test]
    fn sum_chrom_corr_uses_all_allocated_chromosomes() {
        let mut chrom_corr = vec![vec![Corr::default(); 1]; 25];
        chrom_corr[23][0].s0 = 5.0;
        let global = sum_chrom_corr(&chrom_corr, 1);
        assert_eq!(global[0].s0, 5.0);
    }

    #[test]
    fn run_dates_keeps_current_dir_when_helper_fit_fails() {
        let _guard = cwd_lock().lock().unwrap();
        let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/toy");
        let data_dir = tempfile::tempdir().unwrap();
        let work_dir = tempfile::tempdir().unwrap();
        for name in ["toy.geno", "toy.snp", "toy.ind", "poplist.txt"] {
            fs::copy(fixture.join(name), data_dir.path().join(name)).unwrap();
        }
        fs::write(
            data_dir.path().join("par.dates"),
            "genotypename: toy.geno\n\
             snpname: toy.snp\n\
             indivname: toy.ind\n\
             poplistname: poplist.txt\n\
             admixpop: Mix\n\
             output: Toy.out\n\
             binsize: 0.005\n\
             maxdis: 0.025\n\
             seed: 77\n\
             runmode: 1\n\
             checkmap: NO\n\
             numchrom: 2\n\
             qbin: 0\n\
             jackknife: YES\n\
             runfit: YES\n\
             afffit: YES\n\
             lovalfit: 2.5\n",
        )
        .unwrap();

        let original = std::env::current_dir().unwrap();
        std::env::set_current_dir(work_dir.path()).unwrap();
        let result = run_dates(&data_dir.path().join("par.dates"), false);
        let final_dir = std::env::current_dir().unwrap();
        std::env::set_current_dir(original).unwrap();

        let err = result.unwrap_err().to_string();
        assert!(
            err.contains("at least two usable rows"),
            "unexpected error: {err}"
        );
        assert_eq!(final_dir, work_dir.path());
    }
}
