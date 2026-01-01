use hypercube_nets::count_hypercube_nets;
use std::time::Instant;

use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() {
    let max_n = 48; // Adjust this number as needed

    println!(
        "{:<10} | {:<15} | {:<30}",
        "Dimension", "Time", "Number of Nets"
    );
    println!("{:-<10}-|-{:-<15}-|-{:-<30}", "", "", "");

    for n in 1..=max_n {
        let start = Instant::now();
        let nets = count_hypercube_nets(n);
        let duration = start.elapsed();

        println!("{:<10} | {:<15?} | {}", n, duration, nets);
    }
}
