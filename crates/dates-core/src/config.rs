//! Legacy DATES parameter parsing.
//!
//! The original C runtime reads `key: value` parameter files and then performs a
//! light string-substitution pass (`dostrsub`) before extracting typed values.
//! This module preserves that shape while replacing the global `getpars`
//! interface with an owned, testable representation.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow, bail};

/// Parsed legacy parameter file with resolved values.
#[derive(Clone, Debug)]
pub struct LegacyParamFile {
    source: PathBuf,
    entries: BTreeMap<String, String>,
}

/// Typed configuration for the `dates` runtime.
#[derive(Clone, Debug)]
pub struct DatesParams {
    pub par_path: PathBuf,
    pub genotypename: String,
    pub snpname: String,
    pub indivname: String,
    pub poplistname: Option<String>,
    pub weightname: Option<String>,
    pub timeoffsetname: Option<String>,
    pub admixpop: Option<String>,
    pub admixlist: Option<String>,
    pub badsnpname: Option<String>,
    pub output: Option<String>,
    pub runmode: i32,
    pub ldmode: bool,
    pub seed: u64,
    pub nxlim: i32,
    pub minparentcount: i32,
    pub norun: bool,
    pub flatmode: bool,
    pub zdipcorrmode: bool,
    pub jackknife: bool,
    pub maxdis: f64,
    pub binsize: f64,
    pub chithresh: f64,
    pub checkmap: bool,
    pub qbin: usize,
    pub runfit: bool,
    pub afffit: bool,
    pub lovalfit: f64,
    pub chrom: Option<i32>,
    pub nochrom: Option<i32>,
}

impl LegacyParamFile {
    /// Load and resolve a legacy DATES parameter file.
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let raw = fs::read_to_string(path)
            .with_context(|| format!("failed to read parameter file {}", path.display()))?;
        let mut entries = BTreeMap::new();
        for line in raw.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let Some((key, value)) = trimmed.split_once(':') else {
                continue;
            };
            let key = key.trim().to_owned();
            let value = strip_inline_comment(value).trim().to_owned();
            entries.insert(key, value);
        }
        let resolved = resolve_entries(&entries);
        Ok(Self {
            source: path.to_path_buf(),
            entries: resolved,
        })
    }

    /// Return the source path that produced this parameter file.
    pub fn source(&self) -> &Path {
        &self.source
    }

    /// Return the resolved string value for `key`.
    pub fn get(&self, key: &str) -> Option<&str> {
        self.entries.get(key).map(String::as_str)
    }

    /// Return the resolved string value for `key`, cloned into a new string.
    pub fn get_string(&self, key: &str) -> Option<String> {
        self.get(key).map(ToOwned::to_owned)
    }

    /// Return the resolved boolean value for `key`, accepting legacy YES/NO.
    pub fn get_bool(&self, key: &str) -> Result<Option<bool>> {
        match self.get(key) {
            None => Ok(None),
            Some(value) => Ok(Some(parse_legacy_bool(value)?)),
        }
    }

    /// Return the resolved integer value for `key`.
    pub fn get_i32(&self, key: &str) -> Result<Option<i32>> {
        match self.get(key) {
            None => Ok(None),
            Some(value) => {
                Ok(Some(value.parse::<i32>().with_context(|| {
                    format!("invalid integer for key {key}: {value}")
                })?))
            }
        }
    }

    /// Return the resolved floating-point value for `key`.
    pub fn get_f64(&self, key: &str) -> Result<Option<f64>> {
        match self.get(key) {
            None => Ok(None),
            Some(value) => {
                Ok(Some(value.parse::<f64>().with_context(|| {
                    format!("invalid float for key {key}: {value}")
                })?))
            }
        }
    }

    /// Return a copy of the internal resolved map.
    pub fn as_map(&self) -> &BTreeMap<String, String> {
        &self.entries
    }
}

impl DatesParams {
    /// Parse a `dates` parameter file into the typed runtime configuration.
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let params = LegacyParamFile::load(path.as_ref())?;
        let genotypename = required(&params, "genotypename")?;
        let snpname = required(&params, "snpname")?;
        let indivname = required(&params, "indivname")?;
        Ok(Self {
            par_path: params.source().to_path_buf(),
            genotypename,
            snpname,
            indivname,
            poplistname: params.get_string("poplistname"),
            weightname: params.get_string("weightname"),
            timeoffsetname: params
                .get_string("timeoffsetname")
                .or_else(|| params.get_string("timeoffset")),
            admixpop: params.get_string("admixpop"),
            admixlist: params.get_string("admixlist"),
            badsnpname: params.get_string("badsnpname"),
            output: params.get_string("output"),
            runmode: params.get_i32("runmode")?.unwrap_or(0),
            ldmode: params.get_bool("ldmode")?.unwrap_or(false),
            seed: params.get_i32("seed")?.unwrap_or(0).max(0) as u64,
            nxlim: params.get_i32("nxlim")?.unwrap_or(100000),
            minparentcount: params.get_i32("minparentcount")?.unwrap_or(10),
            norun: params.get_bool("norun")?.unwrap_or(false),
            flatmode: params.get_bool("flatmode")?.unwrap_or(false),
            zdipcorrmode: params.get_bool("zdipcorrmode")?.unwrap_or(false),
            jackknife: params.get_bool("jackknife")?.unwrap_or(true),
            maxdis: params.get_f64("maxdis")?.unwrap_or(0.05),
            binsize: params.get_f64("binsize")?.unwrap_or(0.0005),
            chithresh: params.get_f64("chithresh")?.unwrap_or(-1.0),
            checkmap: params.get_bool("checkmap")?.unwrap_or(true),
            qbin: params.get_i32("qbin")?.unwrap_or(0).max(0) as usize,
            runfit: params.get_bool("runfit")?.unwrap_or(true),
            afffit: params.get_bool("afffit")?.unwrap_or(true),
            lovalfit: params.get_f64("lovalfit")?.unwrap_or(0.45),
            chrom: params.get_i32("chrom")?,
            nochrom: params.get_i32("nochrom")?,
        })
    }
}

/// Parse legacy YES/NO or truthy/falsy values.
pub fn parse_legacy_bool(value: &str) -> Result<bool> {
    match value.trim().to_ascii_uppercase().as_str() {
        "YES" | "TRUE" | "1" => Ok(true),
        "NO" | "FALSE" | "0" => Ok(false),
        other => bail!("invalid legacy boolean: {other}"),
    }
}

/// Infer the output prefix used by the Perl helper wrappers.
pub fn infer_output_prefix(
    params: &LegacyParamFile,
    admix_override: Option<&str>,
) -> Result<String> {
    if let Some(output) = params.get("output") {
        return Ok(strip_suffix(output, ".out").to_owned());
    }
    let name = admix_override
        .map(ToOwned::to_owned)
        .or_else(|| params.get_string("admixpop"))
        .ok_or_else(|| anyhow!("could not infer output prefix from output: or admixpop:"))?;
    Ok(name)
}

fn required(params: &LegacyParamFile, key: &str) -> Result<String> {
    params
        .get_string(key)
        .ok_or_else(|| anyhow!("missing required key {key}:"))
}

fn resolve_entries(entries: &BTreeMap<String, String>) -> BTreeMap<String, String> {
    let mut resolved = entries.clone();
    for _ in 0..8 {
        let snapshot = resolved.clone();
        let mut changed = false;
        let mut keys: Vec<_> = snapshot.keys().cloned().collect();
        keys.sort_by_key(|value| usize::MAX - value.len());
        for (name, value) in &snapshot {
            let mut next = value.clone();
            for key in &keys {
                if key == name {
                    continue;
                }
                let replacement = snapshot.get(key).unwrap();
                next = next.replace(&format!("${{{key}}}"), replacement);
                next = replace_legacy_token(&next, key, replacement);
            }
            changed |= next != *value;
            resolved.insert(name.clone(), next);
        }
        if !changed {
            break;
        }
    }
    resolved
}

fn replace_legacy_token(value: &str, key: &str, replacement: &str) -> String {
    let mut out = String::new();
    let bytes = value.as_bytes();
    let key_bytes = key.as_bytes();
    let mut index = 0;
    while index < bytes.len() {
        if bytes[index..].starts_with(key_bytes) {
            let prev_ok = index == 0
                || !bytes[index - 1].is_ascii_alphanumeric()
                    && bytes[index - 1] != b'_'
                    && bytes[index - 1] != b'/';
            let next = index + key_bytes.len();
            let next_ok = next == bytes.len()
                || bytes[next] == b'/'
                || !bytes[next].is_ascii_alphanumeric() && bytes[next] != b'_';
            if prev_ok && next_ok {
                out.push_str(replacement);
                index = next;
                continue;
            }
        }
        out.push(bytes[index] as char);
        index += 1;
    }
    out
}

fn strip_inline_comment(value: &str) -> &str {
    value.split('#').next().unwrap_or(value)
}

fn strip_suffix<'a>(value: &'a str, suffix: &str) -> &'a str {
    value.strip_suffix(suffix).unwrap_or(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_legacy_bool_values() {
        assert!(parse_legacy_bool("YES").unwrap());
        assert!(!parse_legacy_bool("no").unwrap());
    }

    #[test]
    fn resolves_dir_substitution() {
        let mut raw = BTreeMap::new();
        raw.insert("DIR".to_owned(), "./data".to_owned());
        raw.insert("indivname".to_owned(), "DIR/family.ind".to_owned());
        let resolved = resolve_entries(&raw);
        assert_eq!(resolved["indivname"], "./data/family.ind");
    }
}
