use std::io::Write;
use std::process::exit;
use std::time::{Duration, Instant};

use system::matrix::{enumerable, matrix_count, MAX_CHANNELS};
use system::multiverse::{exhaustive, sample, sorted, Execution, Exhaustive, Results, Sampled};
use system::universe::{Rules, DEFAULT_MAX_STEPS, MAX_STEPS_CAP};

const USAGE: &str = "Usage: system [options]

  --size <n>      elements per universe (default 2; channels = 2n, matrices = 2^(4n^2); max 8)
  --steps <n>     steps before classification (default 10, max 32)
  --samples <n>   simulate n random matrices instead of every matrix
  --seed <n>      PRNG seed for --samples (default 1)
  --offset <n>    exhaustive: first matrix index (shard start)
  --stride <n>    exhaustive: index stride (shard count)
  --threads <n>   worker threads (default: all cores)
  --fixed-cycle   also compare the final step with the previous one (see README)
  --quiet         no progress output
  --json          print results as JSON
  -h, --help      this text";

struct Args {
    size: usize,
    rules: Rules,
    samples: Option<u64>,
    seed: u64,
    offset: u128,
    stride: u128,
    threads: Option<usize>,
    quiet: bool,
    json: bool,
}

fn parse_args() -> Result<Args, String> {
    let mut args = Args {
        size: 2,
        rules: Rules::default(),
        samples: None,
        seed: 1,
        offset: 0,
        stride: 1,
        threads: None,
        quiet: false,
        json: false,
    };
    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut value = |name: &str| -> Result<String, String> {
            it.next().ok_or_else(|| format!("{name} requires a value"))
        };
        match flag.as_str() {
            "-h" | "--help" => {
                println!("{USAGE}");
                exit(0);
            }
            "--size" => args.size = num(&value("--size")?)?,
            "--steps" => args.rules.max_steps = num(&value("--steps")?)?,
            "--samples" => args.samples = Some(num(&value("--samples")?)?),
            "--seed" => args.seed = num(&value("--seed")?)?,
            "--offset" => args.offset = num(&value("--offset")?)?,
            "--stride" => args.stride = num(&value("--stride")?)?,
            "--threads" => args.threads = Some(num(&value("--threads")?)?),
            "--fixed-cycle" => args.rules.skip_previous = false,
            "--quiet" => args.quiet = true,
            "--json" => args.json = true,
            other => return Err(format!("unknown argument {other}")),
        }
    }
    if args.size < 1 || args.size * 2 > MAX_CHANNELS {
        return Err(format!("--size must be between 1 and {}", MAX_CHANNELS / 2));
    }
    if args.rules.max_steps < 1 || args.rules.max_steps > MAX_STEPS_CAP {
        return Err(format!(
            "--steps must be between 1 and {MAX_STEPS_CAP} (default {DEFAULT_MAX_STEPS})"
        ));
    }
    if args.stride < 1 {
        return Err("--stride must be at least 1".into());
    }
    Ok(args)
}

fn num<T: std::str::FromStr>(raw: &str) -> Result<T, String> {
    raw.parse().map_err(|_| format!("invalid number: {raw}"))
}

/// Dispatch a runtime `size` to the const-generic implementation.
macro_rules! with_channels {
    ($size:expr, $c:ident => $body:expr) => {
        match $size {
            1 => {
                const $c: usize = 2;
                $body
            }
            2 => {
                const $c: usize = 4;
                $body
            }
            3 => {
                const $c: usize = 6;
                $body
            }
            4 => {
                const $c: usize = 8;
                $body
            }
            5 => {
                const $c: usize = 10;
                $body
            }
            6 => {
                const $c: usize = 12;
                $body
            }
            7 => {
                const $c: usize = 14;
                $body
            }
            8 => {
                const $c: usize = 16;
                $body
            }
            _ => unreachable!("size validated"),
        }
    };
}

fn main() {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("error: {e}\n\n{USAGE}");
            exit(2);
        }
    };
    let channels = args.size * 2;
    let mut exec = Execution::default();
    if let Some(t) = args.threads {
        exec.threads = t.max(1);
    }
    exec.progress = if args.quiet {
        None
    } else {
        Some(Duration::from_secs(5))
    };

    let started = Instant::now();
    let (results, label): (Results, String) = if let Some(samples) = args.samples {
        let job = Sampled {
            samples,
            seed: args.seed,
        };
        let r = with_channels!(args.size, C => sample::<C>(job, args.rules, exec));
        (
            r,
            format!("{samples} sampled matrices (seed {})", args.seed),
        )
    } else {
        if !enumerable(channels) {
            eprintln!(
                "error: size {} has 2^{} matrices; exhaustive mode needs at most 2^64. Use --samples.",
                args.size,
                channels * channels
            );
            exit(1);
        }
        let job = Exhaustive {
            offset: args.offset,
            stride: args.stride,
        };
        let r = with_channels!(args.size, C => exhaustive::<C>(job, args.rules, exec));
        let mut label = format!(
            "{} of {} matrices",
            job.count(channels),
            matrix_count(channels)
        );
        if args.stride > 1 || args.offset > 0 {
            label.push_str(&format!(
                " (offset {}, stride {})",
                args.offset, args.stride
            ));
        }
        (r, label)
    };
    let elapsed = started.elapsed();

    let rows = sorted(&results);
    let total: u64 = rows.iter().map(|&(_, n)| n).sum();
    let mut out = String::new();
    if args.json {
        out.push_str("{\n");
        out.push_str(&format!("  \"size\": {},\n", args.size));
        out.push_str(&format!("  \"steps\": {},\n", args.rules.max_steps));
        out.push_str(&format!(
            "  \"skip_previous\": {},\n",
            args.rules.skip_previous
        ));
        out.push_str(&format!("  \"total\": {total},\n"));
        out.push_str(&format!("  \"elapsed_ms\": {},\n", elapsed.as_millis()));
        out.push_str("  \"results\": {\n");
        for (i, (sig, n)) in rows.iter().enumerate() {
            out.push_str(&format!(
                "    \"{sig}\": {n}{}\n",
                if i + 1 < rows.len() { "," } else { "" }
            ));
        }
        out.push_str("  }\n}\n");
    } else {
        out.push_str(&format!(
            "size {}, {} steps, {label}, {total} simulated in {:.3}s on {} threads\n",
            args.size,
            args.rules.max_steps,
            elapsed.as_secs_f64(),
            exec.threads
        ));
        let width = rows.iter().map(|(s, _)| s.len()).max().unwrap_or(0);
        for (sig, n) in &rows {
            out.push_str(&format!(
                "  {sig:<width$}  {n:>12}  {:>6.2}%\n",
                100.0 * *n as f64 / total.max(1) as f64
            ));
        }
    }
    // A closed pipe (e.g. `| head`) is not an error worth a panic.
    let _ = std::io::stdout().write_all(out.as_bytes());
    let _ = std::io::stdout().flush();
}
