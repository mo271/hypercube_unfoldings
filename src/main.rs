use hypercube_nets::count_hypercube_nets;
use std::env;
use std::time::Instant;

use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() {
    let args: Vec<String> = env::args().collect();

    let (start_n, end_n) = match args.len() {
        1 => (1, u32::MAX),
        2 => {
            let n = args[1].parse().expect("Invalid number");
            (n, n)
        }
        3 => {
            let start = args[1].parse().expect("Invalid start number");
            let end = args[2].parse().expect("Invalid end number");
            (start, end)
        }
        _ => {
            eprintln!("Usage: {} [n] OR [start] [end]", args[0]);
            std::process::exit(1);
        }
    };

    println!(
        "{:<10} | {:<15} | {:<30}",
        "Dimension", "Time", "Number of Nets"
    );
    println!("{:-<10}-|-{:-<15}-|-{:-<30}", "", "", "");

    for n in start_n..=end_n {
        let start = Instant::now();
        let nets = count_hypercube_nets(n);
        let duration = start.elapsed();

        println!("{:<10} | {:<15?} | {}", n, duration, nets);
    }
}
