use std::io::Write;
use std::process::exit;
use std::time::{Duration, Instant};

use dynamics_engine::exact::{
    analyze_clamped_traced, analyze_traced, matrix_orbits, reduced_orbits, Budget, Clamped, Fate,
    Reduced,
};
use dynamics_engine::experiments::{
    basin_map, basin_scan, basins, census_clamped, census_exhaustive, census_sampled, class_name,
    BasinMap, Basins, Census,
};
use dynamics_engine::matrix::{enumerable, matrix_count, MAX_CHANNELS};
use dynamics_engine::multiverse::{exhaustive, sample, sorted, Execution, Exhaustive, Results, Sampled};
use dynamics_engine::universe::{Rules, DEFAULT_MAX_STEPS, MAX_STEPS_CAP};

const USAGE: &str = "Dynamics Engine: exact fates of the discrete dynamical systems
v(t+1) = v(t) + d(sign v(t)).  https://mukos.github.io/dynamics-engine/

Usage: dynamics-engine <command> [options]

Commands:
  show <id>    one universe: its drift rows, every sign change and its fate
  exact        fate census of every universe of one size (--samples for N >= 4)
  basins       fates over every start state in a box, per universe
  basinmap     the attractor each start in a box reaches, for one universe (--index)
  basinscan    basin maps and pattern features of every canonical universe
  orbits       symmetry orbit counts (Burnside)
  run          the original 10-step signature histogram (legacy model)

Examples:
  dynamics-engine show 3-27411247        a spiral whose loop doubles each lap
  dynamics-engine show 2-613             the period-10 cycle from the website
  dynamics-engine exact --size 3         all 387,420,489 three-element universes
  dynamics-engine exact --size 4 --samples 100000000 --seed 7

A universe id is N-index, the same as the website's #u= links. --size N --index I also works.

Options:
  --size <n>        elements per universe (default 2)
  --index <i>       show / basinmap: universe index (base-3 drift rows, see exact.rs)
  --start <v,v,..>  show: initial state instead of the all-fire start
  --samples <n>     random universes instead of all (exact, basins, basinscan, run)
  --seed <n>        PRNG seed for --samples (default 1)
  --max-events <n>  sign changes before giving up (default 10000)
  --spiral <n>      growing returns needed to call a spiral, 0 disables (default 8)
  --radius <r>      basins / basinmap / basinscan: start box [-r, r]^n (default 2)
  --model <m>       sign (default) or clamped (create/destroy with positivity)
  --json            machine-readable output (exact, basins, run)
  --heuristic       exact: also compare with the 10-step heuristic (much slower)
  --threads <n>     worker threads (default: all cores)
  --quiet           no progress output
  -V, --version     print the version
  -h, --help        this text

run (legacy) options:
  --steps <n>       steps before classification (default 10, max 32)
  --offset <n>      exhaustive: first matrix index (shard start)
  --stride <n>      exhaustive: index stride (shard count)
  --fixed-cycle     also compare the final step with the previous one (see README)";

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
        heuristic: false,
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
    if std::env::args().len() == 1 {
        println!("{USAGE}");
        exit(0);
    }
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
            "-V" | "--version" => {
                println!("dynamics-engine {}", env!("CARGO_PKG_VERSION"));
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
            "--heuristic" => args.heuristic = true,
            "--no-heuristic" => args.heuristic = false, // the default; kept for old scripts
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
            id if args.command == "show" && !id.starts_with('-') => {
                // "3-27411247" (size-index, as in the website's #u= links) or a bare index.
                match id.split_once('-') {
                    Some((n, i)) => {
                        args.size = num(n)?;
                        args.index = num(i)?;
                    }
                    None => args.index = num(id)?,
                }
            }
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
    if args.model == "sign" && args.size <= 4 && ["show", "basinmap"].contains(&args.command.as_str()) {
        let count = 3u64.pow(2 * (args.size * args.size) as u32);
        if args.index >= count {
            return Err(format!(
                "N = {} has {count} universes, so the index must be below {count}",
                args.size
            ));
        }
    }
    if ![
        "run",
        "exact",
        "basins",
        "orbits",
        "show",
        "basinmap",
        "basinscan",
    ]
    .contains(&args.command.as_str())
    {
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
            eprintln!("error: {e}\nRun dynamics-engine --help for usage.");
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
        "basinmap" => return cmd_basinmap(&args),
        "basinscan" => return cmd_basinscan(&args, exec.threads),
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
        let p = *r as f64 / c.reduced.max(1) as f64;
        let ci = if args.samples.is_some() {
            format!(
                "  ±{:.4}%",
                100.0 * 1.96 * (p * (1.0 - p) / c.reduced.max(1) as f64).sqrt()
            )
        } else {
            String::new()
        };
        out.push_str(&format!(
            "  {kind:<10} {r:>14} {:>7.3}%{ci}   {m:>20} {:>7.3}%\n",
            pct(*r as u128, c.reduced as u128),
            pct(*m, c.matrices)
        ));
    }
    if args.samples.is_some() {
        out.push_str("  (± = 95% confidence half-width of the reduced-universe fraction)\n");
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
    if args.json {
        out = census_json(args, &c);
    }
    let _ = std::io::stdout().write_all(out.as_bytes());
}

/// Machine-readable census: fate totals and every fine class with a representative universe.
fn census_json(args: &Args, c: &Census) -> String {
    let rows = |m: &[Vec<i8>]| {
        m.iter()
            .map(|row| {
                format!(
                    "[{}]",
                    row.iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                )
            })
            .collect::<Vec<_>>()
            .join(",")
    };
    let mut o = String::new();
    o.push_str(&format!(
        "{{\n  \"model\": \"{}\",\n  \"size\": {},\n  \"reduced\": {},\n  \"matrices\": {},\n  \"sampled\": {},\n",
        args.model,
        args.size,
        c.reduced,
        c.matrices,
        args.samples.is_some()
    ));
    o.push_str("  \"fates\": {");
    o.push_str(
        &c.fates
            .iter()
            .map(|(k, (r, m))| format!("\"{k}\": [{r}, {m}]"))
            .collect::<Vec<_>>()
            .join(", "),
    );
    o.push_str("},\n  \"classes\": [\n");
    let n = c.classes.len();
    for (i, (&(kind, param), &(count, index))) in c.classes.iter().enumerate() {
        let (pos, neg): (Vec<Vec<i8>>, Vec<Vec<i8>>) = if args.model == "clamped" {
            with_channels!(args.size, C => {
                let u = Clamped::<{ C / 2 }>::from_index(index);
                (u.rows.iter().map(|r| r.to_vec()).collect(), Vec::new())
            })
        } else {
            with_channels!(args.size, C => {
                let w = Reduced::<{ C / 2 }>::from_index(index);
                (w.pos.iter().map(|r| r.to_vec()).collect(), w.neg.iter().map(|r| r.to_vec()).collect())
            })
        };
        o.push_str(&format!(
            "    {{\"kind\": \"{}\", \"param\": {param}, \"count\": {count}, \"index\": {index}, \"pos\": [{}], \"neg\": [{}]}}{}\n",
            class_name(kind),
            rows(&pos),
            rows(&neg),
            if i + 1 < n { "," } else { "" }
        ));
    }
    o.push_str("  ]\n}\n");
    o
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
    if args.json {
        out = basins_json(args, &r, starts);
    }
    let _ = std::io::stdout().write_all(out.as_bytes());
}

fn basins_json(args: &Args, r: &Basins, starts: u64) -> String {
    let rows = |m: &Vec<Vec<i8>>| {
        m.iter()
            .map(|row| {
                format!(
                    "[{}]",
                    row.iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                )
            })
            .collect::<Vec<_>>()
            .join(",")
    };
    let mut o = String::new();
    o.push_str(&format!(
        "{{\n  \"size\": {},\n  \"radius\": {},\n  \"universes\": {},\n  \"starts_per_universe\": {},\n  \"sampled\": {},\n",
        args.size,
        args.radius,
        r.reduced,
        starts / r.reduced.max(1),
        args.samples.is_some()
    ));
    o.push_str("  \"by_start\": {");
    o.push_str(
        &r.by_start
            .iter()
            .map(|(k, n)| format!("\"{k}\": {n}"))
            .collect::<Vec<_>>()
            .join(", "),
    );
    o.push_str("},\n  \"kind_sets\": [\n");
    let mut sets: Vec<_> = r.kind_sets.iter().collect();
    sets.sort_by(|a, b| a.1.cmp(b.1).then_with(|| a.0.cmp(b.0)));
    for (i, (k, n)) in sets.iter().enumerate() {
        let rep = &r.representatives[*k];
        let starts = rep
            .starts
            .iter()
            .map(|(v, kind, fate)| {
                format!(
                    "{{\"start\": [{}], \"kind\": \"{kind}\", \"fate\": \"{fate}\"}}",
                    v.iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                )
            })
            .collect::<Vec<_>>()
            .join(", ");
        o.push_str(&format!(
            "    {{\"kinds\": \"{k}\", \"count\": {n}, \"index\": {}, \"pos\": [{}], \"neg\": [{}], \"starts\": [{starts}]}}{}\n",
            rep.index,
            rows(&rep.pos),
            rows(&rep.neg),
            if i + 1 < sets.len() { "," } else { "" }
        ));
    }
    o.push_str("  ]\n}\n");
    o
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
            .join(" ")
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
            out.push_str(&format!("fate: {}\n", fate_text(&fate)));
            if args.model == "sign" && args.start.is_none() && (2..=4).contains(&args.size) {
                out.push_str(&format!(
                    "watch it: https://mukos.github.io/dynamics-engine/#u={}-{}\n",
                    args.size, args.index
                ));
            }
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

fn fate_text(fate: &Fate) -> String {
    match *fate {
        Fate::Fixed { step } => format!("fixed, the velocity reaches zero at step {step}"),
        Fate::Cycle { transient, period } => {
            format!("cycle of period {period}, repeating from step {transient}")
        }
        Fate::Ray { step } => format!("ray, escaping in a straight line from step {step}"),
        Fate::Helix { step, period, shift_norm } => format!(
            "helix, from step {step} a {period}-step segment repeats shifted by a vector of size {shift_norm}"
        ),
        Fate::Spiral { returns, max_norm } => format!(
            "spiral (empirical), {returns} growing returns in a row, |v| up to {max_norm}"
        ),
        Fate::Undecided { events, max_norm } => {
            format!("undecided, no certificate after {events} sign changes, |v| up to {max_norm}")
        }
    }
}

fn attractor_text(a: &dynamics_engine::experiments::Attractor) -> String {
    let v =
        a.id.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",");
    match a.kind {
        0 => format!("fixed at ({v})"),
        1 => format!("cycle through ({v})"),
        2 => format!("ray with drift ({v})"),
        3 => format!("helix with shift ({v})"),
        4 => "spiral".into(),
        _ => "undecided".into(),
    }
}

fn basin_json(m: &BasinMap, size: usize) -> String {
    let attractors = m
        .attractors
        .iter()
        .map(|a| {
            format!(
                "{{\"kind\": \"{}\", \"id\": [{}]}}",
                class_name(a.kind),
                a.id.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .collect::<Vec<_>>()
        .join(", ");
    let cells = m
        .cells
        .iter()
        .map(|c| c.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let f = &m.features;
    format!(
        "{{\"index\": {}, \"size\": {size}, \"radius\": {}, \"features\": {{\"attractors\": {}, \"moving\": {}, \"sinks\": {}, \"inert\": {}, \"kinds\": {}, \"boundary\": {:.4}, \"conic\": {:.4}, \"stripes\": {:.4}, \"share\": [{}]}}, \"attractors\": [{attractors}], \"cells\": [{cells}]}}",
        m.index,
        m.radius,
        f.attractors,
        f.moving,
        f.sinks,
        f.inert,
        f.kinds,
        f.boundary,
        f.conic,
        f.stripes,
        f.share.iter().map(|x| format!("{x:.4}")).collect::<Vec<_>>().join(",")
    )
}

fn feature_line(m: &BasinMap) -> String {
    format!(
        "  universe {:<8} moving {:>3}  sinks {:>3}  inert {:>3}  boundary {:.3}  conic {:.3}  stripes {:.3}\n",
        m.index, m.features.moving, m.features.sinks, m.features.inert, m.features.boundary, m.features.conic, m.features.stripes
    )
}

fn cmd_basinmap(args: &Args) {
    let b = budget(args);
    let m: BasinMap = with_channels!(args.size, C => basin_map(&Reduced::<{ C / 2 }>::from_index(args.index), args.radius, b));
    let mut out = String::new();
    if args.json {
        out.push_str(&basin_json(&m, args.size));
        out.push('\n');
    } else {
        let f = &m.features;
        out.push_str(&format!(
            "universe {} (size {}), starts in [-{}, {}]^{}: {} attractors, {} kinds, boundary {:.3}, conic {:.3}, stripes {:.3}\n",
            m.index, args.size, m.radius, m.radius, args.size, f.attractors, f.kinds, f.boundary, f.conic, f.stripes
        ));
        let glyphs: Vec<char> = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
            .chars()
            .collect();
        for (i, a) in m.attractors.iter().enumerate() {
            let share = m.cells.iter().filter(|&&c| c as usize == i).count();
            out.push_str(&format!(
                "  {} {:<34} {:>6} cells\n",
                glyphs.get(i).copied().unwrap_or('#'),
                attractor_text(a),
                share
            ));
        }
        if args.size == 2 {
            let side = (2 * m.radius + 1) as usize;
            out.push_str("\nmap (v0 across, v1 down from +radius to -radius):\n");
            for row in (0..side).rev() {
                let line: String = (0..side)
                    .map(|col| {
                        glyphs
                            .get(m.cells[col * side + row] as usize)
                            .copied()
                            .unwrap_or('#')
                    })
                    .collect();
                out.push_str(&format!("  {line}\n"));
            }
        }
    }
    let _ = std::io::stdout().write_all(out.as_bytes());
}

fn cmd_basinscan(args: &Args, threads: usize) {
    let b = budget(args);
    let started = Instant::now();
    let maps: Vec<BasinMap> = with_channels!(args.size, C => basin_scan::<{ C / 2 }>(args.radius, b, threads, args.samples.map(|n| (n, args.seed))));
    let elapsed = started.elapsed().as_secs_f64();
    let mut out = String::new();
    if args.json {
        out.push_str(&format!(
            "{{\n  \"size\": {},\n  \"radius\": {},\n  \"universes\": {},\n  \"maps\": [\n",
            args.size,
            args.radius,
            maps.len()
        ));
        for (i, m) in maps.iter().enumerate() {
            out.push_str("    ");
            out.push_str(&basin_json(m, args.size));
            out.push_str(if i + 1 < maps.len() { ",\n" } else { "\n" });
        }
        out.push_str("  ]\n}\n");
    } else {
        out.push_str(&format!(
            "size {}: {} canonical universes, starts in [-{}, {}]^{}, {elapsed:.1}s\n",
            args.size,
            maps.len(),
            args.radius,
            args.radius,
            args.size
        ));
        let mut v: Vec<&BasinMap> = maps.iter().collect();
        v.sort_by(|a, b| {
            b.features
                .attractors
                .cmp(&a.features.attractors)
                .then(a.index.cmp(&b.index))
        });
        out.push_str("\nmost attractors:\n");
        v.iter()
            .take(8)
            .for_each(|m| out.push_str(&feature_line(m)));
        v.sort_by(|a, b| {
            b.features
                .boundary
                .partial_cmp(&a.features.boundary)
                .unwrap()
                .then(a.index.cmp(&b.index))
        });
        out.push_str("\nmost boundary:\n");
        v.iter()
            .take(8)
            .for_each(|m| out.push_str(&feature_line(m)));
        let mut multi: Vec<&BasinMap> = maps
            .iter()
            .filter(|m| m.features.moving + m.features.sinks > 1)
            .collect();
        multi.sort_by(|a, b| {
            a.features
                .conic
                .partial_cmp(&b.features.conic)
                .unwrap()
                .then(a.index.cmp(&b.index))
        });
        out.push_str("\nleast conic (fate changes along rays from the origin):\n");
        multi
            .iter()
            .take(8)
            .for_each(|m| out.push_str(&feature_line(m)));
        v.sort_by(|a, b| {
            b.features
                .stripes
                .partial_cmp(&a.features.stripes)
                .unwrap()
                .then(a.index.cmp(&b.index))
        });
        out.push_str("\nmost striped:\n");
        v.iter()
            .take(8)
            .for_each(|m| out.push_str(&feature_line(m)));
        let conic_all = maps.iter().filter(|m| m.features.conic >= 0.999).count();
        let single = maps.iter().filter(|m| m.features.attractors == 1).count();
        out.push_str(&format!(
            "\nfully conic: {conic_all} of {}   single attractor: {single}\n",
            maps.len()
        ));
    }
    let _ = std::io::stdout().write_all(out.as_bytes());
}
