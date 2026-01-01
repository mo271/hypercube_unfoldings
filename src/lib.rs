use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};
use num_integer::Integer;

pub fn count_hypercube_nets(n: u32) -> BigUint {
    let mut total_sum: BigUint = Zero::zero();
    let mut p = vec![0u32; (n + 1) as usize];
    let mut m = vec![0u32; (n + 1) as usize];
    let factorials: Vec<BigUint> = (0..=n).map(|i| factorial(i)).collect();
    
    iterate_partitions(n, 1, &mut p, &mut m, &factorials, n, &mut total_sum);
    
    let bn_size = (BigUint::from(2u32).pow(n)) * &factorials[n as usize];
    total_sum / bn_size
}

fn factorial(n: u32) -> BigUint {
    let mut f: BigUint = One::one();
    for i in 1..=n { f *= i; }
    f
}

fn iterate_partitions(
    remaining: u32,
    min_k: u32,
    p: &mut Vec<u32>,
    m: &mut Vec<u32>,
    factorials: &[BigUint],
    n: u32,
    total_sum: &mut BigUint,
) {
    if remaining == 0 {
        process_partition(p, m, factorials, n, total_sum);
        return;
    }

    for k in min_k..=remaining {
        let max_i = remaining / k;
        for i in 1..=max_i {
            for j in 0..=i {
                let current_p = j;
                let current_m = i - j;
                p[k as usize] += current_p;
                m[k as usize] += current_m;
                iterate_partitions(remaining - i * k, k + 1, p, m, factorials, n, total_sum);
                p[k as usize] -= current_p;
                m[k as usize] -= current_m;
            }
        }
    }
}

fn process_partition(
    p: &[u32],
    m: &[u32],
    factorials: &[BigUint],
    n: u32,
    total_sum: &mut BigUint,
) {
    let tr_g = calculate_tr_g(p, m, n);
    
    if tr_g.is_zero() { return; }
    
    let cent_size = calculate_cent_size(p, m, factorials);
    let bn_size = (BigUint::from(2u32).pow(n)) * &factorials[n as usize];
    let cl_size = &bn_size / &cent_size;
    let term = &tr_g * &cl_size;
    
    *total_sum += term;
}

fn calculate_cent_size(p: &[u32], m: &[u32], factorials: &[BigUint]) -> BigUint {
    let mut cent: BigUint = One::one();
    for k in 1..p.len() {
        if p[k] > 0 {
            let base = BigUint::from(2 * k as u32);
            cent *= base.pow(p[k]);
            cent *= &factorials[p[k] as usize];
        }
        if m[k] > 0 {
            // m[k] is negative cycle of length 2k.
            // Centralizer size for (1..k -1..-k) is 2k.
            let base = BigUint::from(2 * k as u32);
            cent *= base.pow(m[k]);
            cent *= &factorials[m[k] as usize];
        }
    }
    cent
}

fn calculate_tr_g(p: &[u32], m: &[u32], n: u32) -> BigUint {
    let a_0: u32;
    if p.len() > 1 && p[1] > 0 {
        a_0 = 1;
    } else if (p.len() > 2 && p[2] > 0) || (m.len() > 1 && m[1] > 0) {
        a_0 = 2;
    } else {
        return Zero::zero();
    }
    
    let t_0 = if a_0 == 1 {
        if p[1] < 2 { return Zero::zero(); }
        let c1 = p[1];
        let term1 = BigUint::from(2 * c1 - 2).pow(c1);
        let term2 = if c1 >= 2 { BigUint::from(2 * c1).pow(c1 - 2) } else { Zero::zero() };
        term1 * term2
    } else {
        let p2 = if p.len() > 2 { p[2] } else { 0 };
        let m1 = if m.len() > 1 { m[1] } else { 0 };
        
        let v = (2 * p2 + m1) as usize;
        if v == 0 { return Zero::zero(); }
        if p2 == 0 { return Zero::zero(); }
        
        let size = v - 1;
        if size == 0 { return Zero::zero(); }
        
        let mut mat = vec![vec![BigInt::zero(); size]; size];
        let num_paired_rows = 2 * p2 as usize;
        
        for r in 0..size {
            for c in 0..size {
                let orig_r = r + 1;
                let orig_c = c + 1;
                
                if r == c {
                    if orig_r < num_paired_rows {
                        mat[r][c] = BigInt::from(2 * v as i64 - 3);
                    } else {
                        mat[r][c] = BigInt::from(2 * v as i64 - 2);
                    }
                } else {
                    let is_pair = if orig_r < num_paired_rows && orig_c < num_paired_rows {
                        (orig_r ^ 1) == orig_c
                    } else {
                        false
                    };
                    if is_pair {
                        mat[r][c] = BigInt::from(-1);
                    } else {
                        mat[r][c] = BigInt::from(-2);
                    }
                }
            }
        }
        
        let det = determinant_bigint(&mut mat);
        if det <= Zero::zero() { return Zero::zero(); }
        BigUint::from(2 * p2) * det.to_biguint().unwrap()
    };
    
    let mut prod: BigUint = One::one();
    
    for a in (a_0 + 1)..=(2 * n) {
        let p_a = if (a as usize) < p.len() { p[a as usize] } else { 0 };
        let m_half = if a % 2 == 0 && ((a/2) as usize) < m.len() { m[(a/2) as usize] } else { 0 };
        let n_a = 2 * p_a + m_half;
        
        if n_a == 0 { continue; }
        
        let mut w_a: BigUint = Zero::zero();
        for b in 1..a {
             if a % b == 0 {
                 let p_b = if (b as usize) < p.len() { p[b as usize] } else { 0 };
                 let m_half_b = if b % 2 == 0 && ((b/2) as usize) < m.len() { m[(b/2) as usize] } else { 0 };
                 let n_b = 2 * p_b + m_half_b;
                 w_a += BigUint::from(b) * BigUint::from(n_b);
             }
        }
        
        let term_add = BigUint::from(a) * BigUint::from(n_a);
        let base1 = &w_a + &term_add;
        let exp1 = n_a - 1 - p_a;
        
        // base2 = base1 - 2 (derived from Matrix Tree Theorem)
        let base2 = if base1 >= BigUint::from(2u32) { &base1 - 2u32 } else { Zero::zero() };
        let exp2 = p_a;
        
        let term = w_a * base1.pow(exp1) * base2.pow(exp2);
        prod *= term;
    }
    
    t_0 * prod
}

fn determinant_bigint(matrix: &mut Vec<Vec<BigInt>>) -> BigInt {
    let n = matrix.len();
    if n == 0 { return One::one(); }
    if n == 1 { return matrix[0][0].clone(); }
    
    for k in 0..n-1 {
        if matrix[k][k].is_zero() {
            let mut swap_row = None;
            for i in k+1..n {
                if !matrix[i][k].is_zero() {
                    swap_row = Some(i);
                    break;
                }
            }
            if let Some(r) = swap_row {
                matrix.swap(k, r);
                for i in 0..n { matrix[k][i] = -matrix[k][i].clone(); }
            } else {
                return Zero::zero();
            }
        }
        for i in k+1..n {
            for j in k+1..n {
                let val = &matrix[i][j] * &matrix[k][k] - &matrix[i][k] * &matrix[k][j];
                let denom = if k == 0 { BigInt::one() } else { matrix[k-1][k-1].clone() };
                if denom.is_zero() { return Zero::zero(); }
                let (q, _) = val.div_rem(&denom);
                matrix[i][j] = q;
            }
        }
    }
    matrix[n-1][n-1].clone()
}
