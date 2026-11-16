//! Small helpers: argument parsing, logging, sizes, JSON strings.

use std::collections::HashMap;
use std::sync::OnceLock;
use std::time::Instant;

pub type Res<T> = Result<T, Box<dyn std::error::Error + Send + Sync>>;

static START: OnceLock<Instant> = OnceLock::new();

pub fn elapsed() -> f64 {
    START.get_or_init(Instant::now).elapsed().as_secs_f64()
}

#[macro_export]
macro_rules! log {
    ($($arg:tt)*) => {
        eprintln!("[cpk {:8.1}s] {}", $crate::util::elapsed(), format!($($arg)*))
    };
}

/// `--key value` or `--key=value` options (repeatable) and positional arguments. Every option
/// takes a value; a bare `-` is positional (stdin).
pub struct Args {
    opts: HashMap<String, Vec<String>>,
    pub pos: Vec<String>,
}

impl Args {
    pub fn parse(raw: &[String]) -> Args {
        let mut opts: HashMap<String, Vec<String>> = HashMap::new();
        let mut pos = vec![];
        let mut i = 0;
        while i < raw.len() {
            let a = &raw[i];
            if a.len() > 2 && a.starts_with("--") {
                let key = a[2..].to_string();
                if let Some((k, v)) = key.split_once('=') {
                    opts.entry(k.to_string()).or_default().push(v.to_string());
                } else if i + 1 < raw.len() && !(raw[i + 1].len() > 2 && raw[i + 1].starts_with("--")) {
                    opts.entry(key).or_default().push(raw[i + 1].clone());
                    i += 1;
                } else {
                    opts.entry(key).or_default().push("true".to_string());
                }
            } else {
                pos.push(a.clone());
            }
            i += 1;
        }
        Args { opts, pos }
    }

    pub fn get(&self, key: &str) -> Option<&str> {
        self.opts.get(key).and_then(|v| v.last()).map(|s| s.as_str())
    }

    pub fn all(&self, key: &str) -> Vec<String> {
        self.opts.get(key).cloned().unwrap_or_default()
    }

    pub fn req(&self, key: &str) -> Res<&str> {
        self.get(key).ok_or_else(|| format!("missing required option --{key}").into())
    }

    pub fn usize_or(&self, key: &str, default: usize) -> Res<usize> {
        match self.get(key) {
            None => Ok(default),
            Some(v) => Ok(parse_size(v)? as usize),
        }
    }

    pub fn f64_or(&self, key: &str, default: f64) -> Res<f64> {
        match self.get(key) {
            None => Ok(default),
            Some(v) => v.parse::<f64>().map_err(|_| format!("--{key}: not a number: {v}").into()),
        }
    }
}

/// Parse an integer with an optional k/m/g suffix (powers of ten): "3.1g" = 3_100_000_000.
pub fn parse_size(s: &str) -> Res<u64> {
    let t = s.trim().to_ascii_lowercase();
    let (num, mult) = match t.chars().last() {
        Some('k') => (&t[..t.len() - 1], 1e3),
        Some('m') => (&t[..t.len() - 1], 1e6),
        Some('g') => (&t[..t.len() - 1], 1e9),
        _ => (&t[..], 1.0),
    };
    let v: f64 = num.parse().map_err(|_| format!("not a size: {s}"))?;
    if v < 0.0 {
        return Err(format!("negative size: {s}").into());
    }
    Ok((v * mult).round() as u64)
}

pub fn threads_from(a: &Args) -> Res<usize> {
    let d = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
    Ok(a.usize_or("threads", d)?.max(1))
}

pub fn json_str(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    o.push('"');
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            '\t' => o.push_str("\\t"),
            c if (c as u32) < 0x20 => o.push_str(&format!("\\u{:04x}", c as u32)),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

/// Peak resident set size of this process in bytes.
pub fn peak_rss_bytes() -> u64 {
    unsafe {
        let mut ru: libc::rusage = std::mem::zeroed();
        if libc::getrusage(libc::RUSAGE_SELF, &mut ru) != 0 {
            return 0;
        }
        let v = ru.ru_maxrss as u64;
        // macOS reports bytes, Linux kilobytes
        if cfg!(target_os = "macos") {
            v
        } else {
            v * 1024
        }
    }
}

/// `label=path` pairs from a repeatable option.
pub fn labelled(vals: &[String], what: &str) -> Res<Vec<(String, String)>> {
    let mut out = vec![];
    for v in vals {
        match v.split_once('=') {
            Some((l, p)) if !l.is_empty() && !p.is_empty() => out.push((l.to_string(), p.to_string())),
            _ => return Err(format!("--{what} expects label=path, got '{v}'").into()),
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizes() {
        assert_eq!(parse_size("31").unwrap(), 31);
        assert_eq!(parse_size("3.1g").unwrap(), 3_100_000_000);
        assert_eq!(parse_size("512M").unwrap(), 512_000_000);
        assert!(parse_size("x").is_err());
    }

    #[test]
    fn args() {
        let raw: Vec<String> = ["--k", "31", "--asm", "a=x.fa", "--asm", "b=y.fa", "--out=o", "r1.fq", "-"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let a = Args::parse(&raw);
        assert_eq!(a.get("k"), Some("31"));
        assert_eq!(a.all("asm"), vec!["a=x.fa".to_string(), "b=y.fa".to_string()]);
        assert_eq!(a.get("out"), Some("o"));
        assert_eq!(a.pos, vec!["r1.fq".to_string(), "-".to_string()]);
    }
}
