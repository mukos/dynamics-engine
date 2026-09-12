use std::io::Write;
use std::process::exit;
use std::time::{Duration, Instant};

use system::exact::{
    analyze_clamped_traced, analyze_traced, matrix_orbits, reduced_orbits, Budget, Clamped, Reduced,
};
use system::experiments::{
    basins, census_clamped, census_exhaustive, census_sampled, Basins, Census,
};
use system::matrix::{enumerable, matrix_count, MAX_CHANNELS};
use system::multiverse::{exhaustive, sample, sorted, Execution, Exhaustive, Results, Sampled};
use system::universe::{Rules, DEFAULT_MAX_STEPS, MAX_STEPS_CAP};

const USAGE: &str =
    "Usage: system [run] [options]      10-step heuristic histogram (original model)
       system exact [options]      exact fates of every reduced universe (3^(2n^2) of them)
       system basins [options]     exact fates over all initial states in a box, per universe
       system orbits [options]     symmetry orbit counts (Burnside)
       system show [options]       print one universe's drift rows and event trajectory

run options:

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
  -h, --help      this text

exact / basins options:
  --size <n>        elements (default 2)
  --samples <n>     random reduced universes instead of all (exact only)
  --seed <n>        PRNG seed for --samples
  --max-events <n>  orthant changes before giving up (default 10000)
  --spiral <n>      growing returns needed to call a spiral, 0 disables (default 8)
  --radius <r>      basins: start box [-r, r]^n (default 2)
  --model <m>       sign (original, default) or clamped (create/destroy with positivity)
  --index <i>       show: reduced-universe index (base-3 encoding, see exact.rs)
  --start <v,v,..>  show: initial state instead of the all-fire start
  --no-heuristic    exact: skip the comparison with the 10-step heuristic
  --threads <n>";

struct Args {
    command: String,
    size: usize,
    max_events: u64,
    spiral: u64,
    radius: i64,
    heuristic: bool,
    model: String,
    index: u64,
    start: Option<Vec<i64>>,
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
        command: "run".into(),
        size: 2,
        max_events: 10_000,
        spiral: 8,
        radius: 2,
        heuristic: true,
        model: "sign".into(),
        index: 0,
        start: None,
        rules: Rules::default(),
        samples: None,
        seed: 1,
        offset: 0,
        stride: 1,
        threads: None,
        quiet: false,
        json: false,
    };
    let mut it = std::env::args().skip(1).peekable();
    if let Some(first) = it.peek() {
        if !first.starts_with('-') {
            args.command = it.next().unwrap();
        }
    }
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
            "--max-events" => args.max_events = num(&value("--max-events")?)?,
            "--spiral" => args.spiral = num(&value("--spiral")?)?,
            "--radius" => args.radius = num(&value("--radius")?)?,
            "--no-heuristic" => args.heuristic = false,
            "--model" => args.model = value("--model")?,
            "--index" => args.index = num(&value("--index")?)?,
            "--start" => {
                let raw = value("--start")?;
                args.start = Some(
                    raw.split(',')
                        .map(|x| num::<i64>(x.trim()))
                        .collect::<Result<_, _>>()?,
                );
            }
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
    if !["sign", "clamped"].contains(&args.model.as_str()) {
        return Err("--model must be sign or clamped".into());
    }
    if !["run", "exact", "basins", "orbits", "show"].contains(&args.command.as_str()) {
        return Err(format!("unknown command {}", args.command));
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

    match args.command.as_str() {
        "exact" => return cmd_exact(&args, exec.threads),
        "basins" => return cmd_basins(&args, exec.threads),
        "orbits" => return cmd_orbits(&args),
        "show" => return cmd_show(&args),
        _ => {}
    }

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

fn budget(args: &Args) -> Budget {
    Budget {
        max_events: args.max_events,
        spiral_returns: args.spiral,
        ..Budget::default()
    }
}

fn pct(n: u128, total: u128) -> f64 {
    100.0 * n as f64 / total.max(1) as f64
}

fn cmd_exact(args: &Args, threads: usize) {
    let started = Instant::now();
    let b = budget(args);
    let c: Census = match (args.model.as_str(), args.samples) {
        ("clamped", _) => with_channels!(args.size, C => census_clamped::<{ C / 2 }>(b, threads)),
        (_, Some(n)) => {
            with_channels!(args.size, C => census_sampled::<{ C / 2 }>(n, args.seed, b, threads, args.heuristic))
        }
        _ => {
            with_channels!(args.size, C => census_exhaustive::<{ C / 2 }>(b, threads, args.heuristic))
        }
    };
    let elapsed = started.elapsed().as_secs_f64();
    let mut out = String::new();
    out.push_str(&format!(
        "model {}, size {}: {} reduced universes ({} bit matrices) in {elapsed:.2}s on {threads} threads\n",
        args.model, args.size, c.reduced, c.matrices
    ));
    out.push_str("\nfate          reduced universes        bit-matrix weighted\n");
    for (kind, (r, m)) in &c.fates {
        out.push_str(&format!(
            "  {kind:<10} {r:>14} {:>6.2}%   {m:>20} {:>6.2}%\n",
            pct(*r as u128, c.reduced as u128),
            pct(*m, c.matrices)
        ));
    }
    if c.heuristic_available() {
        out.push_str("\n10-step heuristic verdict vs exact fate (bit-matrix weighted):\n");
        let kinds = ["fixed", "cycle", "ray", "helix", "spiral", "undecided"];
        out.push_str(&format!("  {:<12}", "heuristic"));
        for k in kinds {
            out.push_str(&format!("{k:>16}"));
        }
        out.push('\n');
        for verdict in ["terminated", "bounded", "growth"] {
            out.push_str(&format!("  {verdict:<12}"));
            for k in kinds {
                let n = c.confusion.get(&(verdict, k)).copied().unwrap_or(0);
                out.push_str(&format!("{:>15.2}%", pct(n, c.matrices)));
            }
            out.push('\n');
        }
    }
    if !c.periods.is_empty() {
        out.push_str("\ncycle periods (reduced universes):\n");
        for (p, n) in &c.periods {
            out.push_str(&format!("  period {p:>4}: {n}\n"));
        }
        out.push_str(&format!(
            "  longest transient+period: {} (reduced index {})\n",
            c.longest_cycle.0, c.longest_cycle.1
        ));
    }
    if !c.certified_at.is_empty() {
        let max_step = c.certified_at.keys().max().copied().unwrap_or(0);
        let within10: u64 = c
            .certified_at
            .iter()
            .filter(|(&s, _)| s <= 10)
            .map(|(_, &n)| n)
            .sum();
        let all: u64 = c.certified_at.values().sum();
        out.push_str(&format!(
            "\nfixed/ray/helix certified: {within10} of {all} by step 10, latest at step {max_step}\n"
        ));
    }
    if !c.spiral_examples.is_empty() {
        out.push_str(&format!(
            "\nspiral examples (reduced indices): {:?}\n",
            c.spiral_examples
        ));
    }
    if !c.undecided_examples.is_empty() {
        out.push_str(&format!(
            "\nundecided examples (reduced indices): {:?}\n",
            c.undecided_examples
        ));
    }
    let _ = std::io::stdout().write_all(out.as_bytes());
}

fn cmd_basins(args: &Args, threads: usize) {
    let started = Instant::now();
    let b = budget(args);
    let r: Basins = with_channels!(args.size, C => basins::<{ C / 2 }>(args.radius, b, threads, args.samples.map(|n| (n, args.seed))));
    let elapsed = started.elapsed().as_secs_f64();
    let mut out = String::new();
    let starts: u64 = r.by_start.values().sum();
    out.push_str(&format!(
        "size {}: {} reduced universes x {} starts in [-{}, {}]^{} in {elapsed:.2}s\n",
        args.size,
        r.reduced,
        starts / r.reduced.max(1),
        args.radius,
        args.radius,
        args.size
    ));
    out.push_str("\nfate over all (universe, start) pairs:\n");
    for (k, n) in &r.by_start {
        out.push_str(&format!(
            "  {k:<10} {n:>12} {:>6.2}%\n",
            pct(*n as u128, starts as u128)
        ));
    }
    out.push_str("\ndistinct fate kinds per universe across starts:\n");
    for (k, n) in &r.kinds_per_universe {
        out.push_str(&format!(
            "  {k} kind(s): {n:>10} {:>6.2}%\n",
            pct(*n as u128, r.reduced as u128)
        ));
    }
    out.push_str("\nfate-kind sets:\n");
    let mut sets: Vec<_> = r.kind_sets.iter().collect();
    sets.sort_by(|a, b| b.1.cmp(a.1));
    for (k, n) in sets {
        out.push_str(&format!(
            "  {k:<40} {n:>10} {:>6.2}%\n",
            pct(*n as u128, r.reduced as u128)
        ));
    }
    let _ = std::io::stdout().write_all(out.as_bytes());
}

fn cmd_orbits(args: &Args) {
    let mut out = String::new();
    out.push_str(
        "size  group   reduced universes   orbits (reduced)   bit matrices   orbits (bit)\n",
    );
    for size in 1..=args.size.min(4) {
        let (ro, mo, rc, mc) = match size {
            1 => (
                reduced_orbits::<1>(),
                matrix_orbits::<1>(),
                Reduced::<1>::count(),
                matrix_count(2),
            ),
            2 => (
                reduced_orbits::<2>(),
                matrix_orbits::<2>(),
                Reduced::<2>::count(),
                matrix_count(4),
            ),
            3 => (
                reduced_orbits::<3>(),
                matrix_orbits::<3>(),
                Reduced::<3>::count(),
                matrix_count(6),
            ),
            _ => (
                reduced_orbits::<4>(),
                matrix_orbits::<4>(),
                Reduced::<4>::count(),
                matrix_count(8),
            ),
        };
        let group = (1u128 << size) * (1..=size as u128).product::<u128>();
        out.push_str(&format!(
            "{size:>4}  {group:>5}   {rc:>17}   {ro:>16}   {mc:>12}   {mo:>12}\n"
        ));
    }
    let _ = std::io::stdout().write_all(out.as_bytes());
}

fn cmd_show(args: &Args) {
    let b = budget(args);
    let mut out = String::new();
    let fmt = |v: &[i64]| {
        v.iter()
            .map(|x| format!("{x:>4}"))
            .collect::<Vec<_>>()
            .join("")
    };
    macro_rules! show {
        ($n:literal) => {{
            let start_of = |default: [i64; $n]| -> [i64; $n] {
                match &args.start {
                    Some(s) if s.len() == $n => std::array::from_fn(|i| s[i]),
                    Some(_) => {
                        eprintln!("--start needs {} values", $n);
                        exit(2)
                    }
                    None => default,
                }
            };
            let mut events = String::new();
            let mut trace = |step: u64, v: &[i64; $n], d: &[i64; $n]| {
                events.push_str(&format!("  step {step:>6}  v = [{} ]  drift = [{} ]\n", fmt(v), fmt(d)));
            };
            let fate = if args.model == "clamped" {
                let c = Clamped::<$n>::from_index(args.index);
                out.push_str(&format!(
                    "clamped universe {} (multiplicity {}), rows (element f fires while v_f > 0):\n",
                    args.index,
                    c.multiplicity()
                ));
                for f in 0..$n {
                    out.push_str(&format!("  row {f}: [{} ]\n", fmt(&c.rows[f].map(|x| x as i64))));
                }
                analyze_clamped_traced(&c, start_of(c.initial_state()), b, &mut trace)
            } else {
                let w = Reduced::<$n>::from_index(args.index);
                out.push_str(&format!(
                    "reduced universe {} (multiplicity {}, canonical: {}):\n",
                    args.index,
                    w.multiplicity(),
                    w.canonical() == w
                ));
                for f in 0..$n {
                    out.push_str(&format!(
                        "  element {f}: pos = [{} ]  neg = [{} ]\n",
                        fmt(&w.pos[f].map(|x| x as i64)),
                        fmt(&w.neg[f].map(|x| x as i64))
                    ));
                }
                analyze_traced(&w, start_of(w.initial_state()), 1, b, &mut trace)
            };
            out.push_str("events:\n");
            out.push_str(&events);
            out.push_str(&format!("fate: {fate:?}\n"));
        }};
    }
    match args.size {
        1 => show!(1),
        2 => show!(2),
        3 => show!(3),
        4 => show!(4),
        _ => {
            eprintln!("show supports sizes 1-4");
            exit(2)
        }
    }
    let _ = std::io::stdout().write_all(out.as_bytes());
}
