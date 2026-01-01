use hypercube_nets::count_hypercube_nets;
use std::time::Instant;

fn main() {
    let max_n = 48; // Adjust this number as needed
    
    println!("{:<10} | {:<30} | {:<15}", "Dimension", "Number of Nets", "Time");
    println!("{:-<10}-|-{:-<30}-|-{:-<15}", "", "", "");
    
    for n in 1..=max_n {
        let start = Instant::now();
        let nets = count_hypercube_nets(n);
        let duration = start.elapsed();
        
        println!("{:<10} | {:<30} | {:?}", n, nets, duration);
    }
}
