//! Plot artifact generation for `dates_plot`.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result};

use crate::dataset::FitRow;

/// Plot description shared by the `.xtxt`, `.ps`, and `.pdf` outputs.
#[derive(Clone, Debug)]
pub struct PlotSpec {
    pub title: String,
    pub x_label: String,
    pub y_label: String,
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
    pub data: Vec<(f64, f64)>,
    pub fit: Vec<(f64, f64)>,
}

impl PlotSpec {
    /// Build a plot spec from fitted rows.
    pub fn from_fit(title: impl Into<String>, rows: &[FitRow], x_min: f64, x_max: f64) -> Self {
        let title = title.into();
        let mut y_min = f64::INFINITY;
        let mut y_max = f64::NEG_INFINITY;
        let data = rows
            .iter()
            .map(|row| (row.distance_cm, row.observed))
            .collect::<Vec<_>>();
        let fit = rows
            .iter()
            .map(|row| (row.distance_cm, row.fitted))
            .collect::<Vec<_>>();
        for (_, value) in data.iter().chain(fit.iter()) {
            y_min = y_min.min(*value);
            y_max = y_max.max(*value);
        }
        let pad = ((y_max - y_min).abs() * 0.1).max(1.0e-6);
        Self {
            title,
            x_label: "Genetic Distance (cM)".to_owned(),
            y_label: "Weighted Covariance".to_owned(),
            x_min,
            x_max,
            y_min: y_min - pad,
            y_max: y_max + pad,
            data,
            fit,
        }
    }
}

/// Write the legacy gnuplot-style `.xtxt` file.
pub fn write_xtxt(path: &Path, spec: &PlotSpec, fit_path: &Path) -> Result<()> {
    let fit_name = fit_path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("fit.out");
    let body = format!(
        "set terminal postscript color\n\
         set title  \"{}\"\n\
         set key top right\n\
         set xlabel  \"{}\"\n\
         set ylabel  \"{}\"\n\
         set xrange [{}:{}]\n\
         plot \"{}\" using 1:2  title \"data\", \"{}\" using 1:3 title \"fit\" with lines\n",
        spec.title, spec.x_label, spec.y_label, spec.x_min, spec.x_max, fit_name, fit_name
    );
    fs::write(path, body).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

/// Write a simple PostScript plot artifact.
pub fn write_ps(path: &Path, spec: &PlotSpec) -> Result<()> {
    let mut out = String::new();
    out.push_str("%!PS-Adobe-3.0\n");
    out.push_str("%%Creator: DATES-rs\n");
    out.push_str("/Helvetica findfont 12 scalefont setfont\n");
    out.push_str(&format!(
        "80 760 moveto ({}) show\n",
        ps_escape(&spec.title)
    ));
    out.push_str("newpath 70 100 moveto 70 500 lineto 550 500 lineto stroke\n");
    out.push_str("newpath 70 100 moveto 550 100 lineto stroke\n");
    out.push_str(&render_ps_series(spec, &spec.data, "0 0 0"));
    out.push_str(&render_ps_series(spec, &spec.fit, "1 0 0"));
    fs::write(path, out).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

/// Write a small vector PDF plot artifact.
pub fn write_pdf(path: &Path, spec: &PlotSpec) -> Result<()> {
    let content = render_pdf_content(spec);
    let objects = vec![
        "1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n".to_owned(),
        "2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n".to_owned(),
        "3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>\nendobj\n".to_owned(),
        "4 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>\nendobj\n".to_owned(),
        format!(
            "5 0 obj\n<< /Length {} >>\nstream\n{}endstream\nendobj\n",
            content.len(),
            content
        ),
    ];
    let mut pdf = String::from("%PDF-1.4\n");
    let mut offsets = vec![0usize];
    for object in &objects {
        offsets.push(pdf.len());
        pdf.push_str(object);
    }
    let xref_offset = pdf.len();
    pdf.push_str(&format!("xref\n0 {}\n", offsets.len()));
    pdf.push_str("0000000000 65535 f \n");
    for offset in offsets.iter().skip(1) {
        pdf.push_str(&format!("{offset:010} 00000 n \n"));
    }
    pdf.push_str(&format!(
        "trailer\n<< /Size {} /Root 1 0 R >>\nstartxref\n{}\n%%EOF\n",
        offsets.len(),
        xref_offset
    ));
    fs::write(path, pdf).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn ps_escape(input: &str) -> String {
    input
        .replace('\\', "\\\\")
        .replace('(', "\\(")
        .replace(')', "\\)")
}

fn render_ps_series(spec: &PlotSpec, series: &[(f64, f64)], color: &str) -> String {
    let mut out = String::new();
    out.push_str(&format!("{color} setrgbcolor\n"));
    if let Some((first_x, first_y)) = series.first() {
        let (x, y) = map_point(spec, *first_x, *first_y);
        out.push_str(&format!("newpath {x:.2} {y:.2} moveto\n"));
        for (distance, value) in series.iter().skip(1) {
            let (x, y) = map_point(spec, *distance, *value);
            out.push_str(&format!("{x:.2} {y:.2} lineto\n"));
        }
        out.push_str("stroke\n");
    }
    out
}

fn render_pdf_content(spec: &PlotSpec) -> String {
    let mut out = String::new();
    out.push_str("BT /F1 14 Tf 80 760 Td ");
    out.push_str(&format!("({}) Tj ET\n", ps_escape(&spec.title)));
    out.push_str("0 0 0 RG 1 w\n");
    out.push_str("70 100 m 70 500 l 550 500 l S\n");
    out.push_str("70 100 m 550 100 l S\n");
    out.push_str(&render_pdf_series(spec, &spec.data, "0 0 0"));
    out.push_str(&render_pdf_series(spec, &spec.fit, "1 0 0"));
    out
}

fn render_pdf_series(spec: &PlotSpec, series: &[(f64, f64)], color: &str) -> String {
    let mut out = String::new();
    out.push_str(&format!("{color} RG 1.5 w\n"));
    if let Some((first_x, first_y)) = series.first() {
        let (x, y) = map_point(spec, *first_x, *first_y);
        out.push_str(&format!("{x:.2} {y:.2} m\n"));
        for (distance, value) in series.iter().skip(1) {
            let (x, y) = map_point(spec, *distance, *value);
            out.push_str(&format!("{x:.2} {y:.2} l\n"));
        }
        out.push_str("S\n");
    }
    out
}

fn map_point(spec: &PlotSpec, x: f64, y: f64) -> (f64, f64) {
    let x_scale = 480.0 / (spec.x_max - spec.x_min);
    let y_scale = 400.0 / (spec.y_max - spec.y_min);
    (
        70.0 + (x - spec.x_min) * x_scale,
        100.0 + (y - spec.y_min) * y_scale,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plot_files_are_written() {
        let dir = tempfile::tempdir().unwrap();
        let rows = vec![
            FitRow {
                distance_cm: 0.5,
                observed: 0.02,
                fitted: 0.018,
                residual: 0.002,
            },
            FitRow {
                distance_cm: 0.6,
                observed: 0.018,
                fitted: 0.017,
                residual: 0.001,
            },
        ];
        let spec = PlotSpec::from_fit("DATES: test", &rows, 0.5, 20.0);
        write_ps(&dir.path().join("plot.ps"), &spec).unwrap();
        write_pdf(&dir.path().join("plot.pdf"), &spec).unwrap();
        write_xtxt(&dir.path().join("plot.xtxt"), &spec, Path::new("fit.out")).unwrap();
    }

    #[test]
    fn xtxt_uses_requested_x_range() {
        let dir = tempfile::tempdir().unwrap();
        let rows = vec![
            FitRow {
                distance_cm: 1.0,
                observed: 0.02,
                fitted: 0.018,
                residual: 0.002,
            },
            FitRow {
                distance_cm: 5.0,
                observed: 0.018,
                fitted: 0.017,
                residual: 0.001,
            },
        ];
        let spec = PlotSpec::from_fit("DATES: test", &rows, 1.0, 5.0);
        let path = dir.path().join("plot.xtxt");
        write_xtxt(&path, &spec, Path::new("fit.out")).unwrap();
        let body = fs::read_to_string(path).unwrap();
        assert!(body.contains("set xrange [1:5]"));
    }

    #[test]
    fn ps_and_pdf_coordinates_follow_requested_x_range() {
        let rows = vec![
            FitRow {
                distance_cm: 1.0,
                observed: 0.02,
                fitted: 0.018,
                residual: 0.002,
            },
            FitRow {
                distance_cm: 5.0,
                observed: 0.018,
                fitted: 0.017,
                residual: 0.001,
            },
        ];
        let default = PlotSpec::from_fit("DATES: test", &rows, 0.5, 20.0);
        let narrowed = PlotSpec::from_fit("DATES: test", &rows, 1.0, 5.0);
        let default_ps = render_ps_series(&default, &default.data, "0 0 0");
        let narrowed_ps = render_ps_series(&narrowed, &narrowed.data, "0 0 0");
        let default_pdf = render_pdf_series(&default, &default.data, "0 0 0");
        let narrowed_pdf = render_pdf_series(&narrowed, &narrowed.data, "0 0 0");

        assert_ne!(default_ps, narrowed_ps);
        assert_ne!(default_pdf, narrowed_pdf);
        assert!(narrowed_ps.contains("70.00"));
        assert!(narrowed_ps.contains("550.00"));
        assert!(narrowed_pdf.contains("70.00"));
        assert!(narrowed_pdf.contains("550.00"));
    }
}
