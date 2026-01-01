use rayon::prelude::*;
use rug::{
    Assign, Integer,
    ops::{Pow, SubFrom},
};

pub fn count_hypercube_nets(n: u32) -> Integer {
    let factorials: Vec<Integer> = (0..=n).map(|i| factorial(i)).collect();
    let bn_size = Integer::from(2u32).pow(n) * &factorials[n as usize];

    // Precompute partitions for all k in 0..=n
    // This is fast enough to be serial or parallel, but generating is recursive.
    // We can just generate them.
    let partitions = generate_all_partitions(n);

    // Precompute data for P partitions in parallel
    // (0..=n).into_par_iter() preserves order in collect.
    let precomputed_p: Vec<Vec<PreP>> = (0..=n)
        .into_par_iter()
        .map(|n_p| {
            partitions[n_p as usize]
                .iter()
                .map(|p_vec| precompute_p(p_vec, &factorials, n))
                .collect()
        })
        .collect();

    // Precompute data for M partitions in parallel
    let precomputed_m: Vec<Vec<PreM>> = (0..=n)
        .into_par_iter()
        .map(|n_m| {
            partitions[n_m as usize]
                .iter()
                .map(|m_vec| precompute_m(m_vec, &factorials, n))
                .collect()
        })
        .collect();

    // Summing up terms
    // We iterate n_p sequentially (or parallel? n is small, 30-100).
    // Parallelizing the inner loops is more important.
    // But we can flatten the work if we want.
    // Strategy: For each n_p, parallelize over p_list.
    // This gives plenty of tasks (sum p(k) over k).

    let total_sum: Integer = (0..=n)
        .into_par_iter()
        .map(|n_p| {
            let n_m = n - n_p;
            let p_list = &precomputed_p[n_p as usize];
            let m_list = &precomputed_m[n_m as usize];

            // Parallelize over p_list
            // Each p_data processes m_list (serial).
            // This is good granularity.
            p_list
                .par_iter()
                .map(|p_data| {
                    let mut local_sum = Integer::ZERO;
                    for m_data in m_list {
                        local_sum += calculate_term(p_data, m_data, n, &bn_size);
                    }
                    local_sum
                })
                .reduce(|| Integer::ZERO, |a, b| a + b)
        })
        .reduce(|| Integer::ZERO, |a, b| a + b);

    total_sum / bn_size
}

struct PreP {
    p: Vec<u32>,
    cent: Integer,
    n_contrib: Vec<u32>,
    w_contrib: Vec<Integer>,
    t0_partial_a0_1: Integer,
}

struct PreM {
    m1: u32,
    cent: Integer,
    n_contrib: Vec<u32>,
    w_contrib: Vec<Integer>,
}

fn generate_all_partitions(n: u32) -> Vec<Vec<Vec<u32>>> {
    let mut result = vec![Vec::new(); (n + 1) as usize];
    result[0].push(vec![0; (n + 1) as usize]);

    for target in 1..=n {
        let mut parts = Vec::new();
        generate_partitions_recursive(target, 1, &mut vec![0; (n + 1) as usize], &mut parts);
        result[target as usize] = parts;
    }
    result
}

fn generate_partitions_recursive(
    remaining: u32,
    min_k: u32,
    current: &mut Vec<u32>,
    result: &mut Vec<Vec<u32>>,
) {
    if remaining == 0 {
        result.push(current.clone());
        return;
    }
    for k in min_k..=remaining {
        let max_i = remaining / k;
        for i in 1..=max_i {
            current[k as usize] += i;
            generate_partitions_recursive(remaining - i * k, k + 1, current, result);
            current[k as usize] -= i;
        }
    }
}

fn precompute_p(p: &[u32], factorials: &[Integer], n: u32) -> PreP {
    let mut cent = Integer::ONE.clone();
    for k in 1..p.len() {
        if p[k] > 0 {
            let base = Integer::from(2 * k as u32);
            cent *= base.pow(p[k]);
            cent *= &factorials[p[k] as usize];
        }
    }

    let mut n_contrib = vec![0u32; (2 * n + 1) as usize];
    let mut w_contrib = vec![Integer::ZERO; (2 * n + 1) as usize];

    for k in 1..p.len() {
        if p[k] > 0 {
            if (k as usize) < n_contrib.len() {
                n_contrib[k] = n_contrib[k] + 2 * p[k];
            }
        }
    }

    for a in 1..w_contrib.len() {
        let mut w = Integer::ZERO.clone();
        for b in 1..a {
            if a % b == 0 {
                if b < p.len() && p[b] > 0 {
                    w += Integer::from(b as u32) * Integer::from(2 * p[b]);
                }
            }
        }
        w_contrib[a] = w;
    }

    let t0_partial_a0_1 = if p.len() > 1 && p[1] >= 2 {
        let c1 = p[1];
        Integer::from(2 * c1 - 2).pow(c1) * Integer::from(2 * c1).pow(c1 - 2)
    } else {
        Integer::ZERO
    };

    PreP {
        p: p.to_vec(),
        cent,
        n_contrib,
        w_contrib,
        t0_partial_a0_1,
    }
}

fn precompute_m(m: &[u32], factorials: &[Integer], n: u32) -> PreM {
    let mut cent = Integer::ONE.clone();
    for k in 1..m.len() {
        if m[k] > 0 {
            let base = Integer::from(2 * k as u32);
            cent *= base.pow(m[k]);
            cent *= &factorials[m[k] as usize];
        }
    }

    let mut n_contrib = vec![0u32; (2 * n + 1) as usize];
    let mut w_contrib = vec![Integer::ZERO; (2 * n + 1) as usize];

    for k in 1..m.len() {
        if m[k] > 0 {
            let a = 2 * k;
            if (a as usize) < n_contrib.len() {
                n_contrib[a] += m[k];
            }
        }
    }

    for a in 1..w_contrib.len() {
        let mut w = Integer::ZERO;
        for b in 1..a {
            if a % b == 0 {
                if b % 2 == 0 {
                    let k = b / 2;
                    if k < m.len() && m[k] > 0 {
                        w += Integer::from(b as u32) * Integer::from(m[k]);
                    }
                }
            }
        }
        w_contrib[a] = w;
    }

    let m1 = if m.len() > 1 { m[1] } else { 0 };

    PreM {
        m1,
        cent,
        n_contrib,
        w_contrib,
    }
}

fn calculate_term(p_data: &PreP, m_data: &PreM, n: u32, bn_size: &Integer) -> Integer {
    let a_0: u32;
    if p_data.p.len() > 1 && p_data.p[1] > 0 {
        a_0 = 1;
    } else if (p_data.p.len() > 2 && p_data.p[2] > 0) || m_data.m1 > 0 {
        a_0 = 2;
    } else {
        return Integer::ZERO;
    }

    let t_0 = if a_0 == 1 {
        p_data.t0_partial_a0_1.clone()
    } else {
        let p2 = if p_data.p.len() > 2 { p_data.p[2] } else { 0 };
        let m1 = m_data.m1;
        let v = (2 * p2 + m1) as usize;

        if v == 0 || p2 == 0 {
            return Integer::ZERO;
        }

        let size = v - 1;
        let mut mat = vec![vec![Integer::ZERO; size]; size];
        let num_paired_rows = 2 * p2 as usize;

        for r in 0..size {
            for c in 0..size {
                let orig_r = r + 1;
                let orig_c = c + 1;
                if r == c {
                    if orig_r < num_paired_rows {
                        mat[r][c] = Integer::from(2 * v as i64 - 3);
                    } else {
                        mat[r][c] = Integer::from(2 * v as i64 - 2);
                    }
                } else {
                    let is_pair = if orig_r < num_paired_rows && orig_c < num_paired_rows {
                        (orig_r ^ 1) == orig_c
                    } else {
                        false
                    };
                    if is_pair {
                        mat[r][c] = Integer::from(-1);
                    } else {
                        mat[r][c] = Integer::from(-2);
                    }
                }
            }
        }
        let det = determinant_bigint(&mut mat);
        if det <= 0 {
            return Integer::ZERO;
        }
        Integer::from(2 * p2) * det
    };

    if t_0.is_zero() {
        return Integer::ZERO;
    }

    let mut prod = Integer::ONE.clone();

    for a in (a_0 + 1)..=(2 * n) {
        let a_idx = a as usize;
        let n_a = if a_idx < p_data.n_contrib.len() {
            p_data.n_contrib[a_idx]
        } else {
            0
        } + if a_idx < m_data.n_contrib.len() {
            m_data.n_contrib[a_idx]
        } else {
            0
        };

        if n_a == 0 {
            continue;
        }

        let w_a = p_data.w_contrib[a_idx].clone() + &m_data.w_contrib[a_idx];

        let p_a = if a_idx < p_data.p.len() {
            p_data.p[a_idx]
        } else {
            0
        };

        let term_add = Integer::from(a) * Integer::from(n_a);
        let base1 = w_a.clone() + &term_add;
        let exp1 = n_a - 1 - p_a;

        let base2 = if base1 >= Integer::from(2u32) {
            base1.clone() - 2u32
        } else {
            Integer::ZERO
        };
        let exp2 = p_a;

        let term = w_a * base1.pow(exp1) * base2.pow(exp2);
        prod *= term;
    }

    let tr_g = t_0 * prod;
    let cent_size = p_data.cent.clone() * &m_data.cent;
    (tr_g * bn_size) / cent_size
}

fn determinant_bigint(matrix: &mut Vec<Vec<Integer>>) -> Integer {
    let n = matrix.len();
    if n == 0 {
        return Integer::ONE.clone();
    }
    if n == 1 {
        return matrix[0][0].clone();
    }

    // Re-use storage to avoid allocations in hot loop.
    let mut val_tmp = Integer::new();
    let mut val = Integer::new();

    for k in 0..n - 1 {
        if matrix[k][k].is_zero() {
            let mut swap_row = None;
            for i in k + 1..n {
                if !matrix[i][k].is_zero() {
                    swap_row = Some(i);
                    break;
                }
            }
            if let Some(r) = swap_row {
                matrix.swap(k, r);
                for i in 0..n {
                    matrix[k][i].sub_from(Integer::ZERO);
                }
            } else {
                return Integer::ZERO;
            }
        }
        for i in k + 1..n {
            for j in k + 1..n {
                val_tmp.assign(&matrix[i][j] * &matrix[k][k]);
                val.assign(&matrix[i][k] * &matrix[k][j]);
                // val = val_tmp - val
                val.sub_from(&val_tmp);
                let (before, after) = matrix.split_at_mut(k);
                let denom = if k == 0 {
                    Integer::ONE
                } else {
                    &before[k - 1][k - 1]
                };
                if denom.is_zero() {
                    return Integer::ZERO;
                }
                after[i - k][j].assign(val.div_exact_ref(denom));
            }
        }
    }
    matrix[n - 1][n - 1].clone()
}

fn factorial(n: u32) -> Integer {
    let mut f = Integer::ONE.clone();
    for i in 1..=n {
        f *= i;
    }
    f
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn test_values() {
        let expected = vec![
            (1, "0"),
            (2, "1"),
            (3, "11"),
            (4, "261"),
            (5, "9694"),
            (6, "502110"),
            (7, "33064966"),
            (8, "2642657228"),
            (9, "248639631948"),
            (10, "26941775019280"),
            (11, "3306075027570423"),
            (12, "453373928307505005"),
            (13, "68734915059053558299"),
            (14, "11418459384326497964902"),
            (15, "2062999819948725194529075"),
            (16, "402798929430911987111828116"),
            (17, "84526877217018050866911342594"),
            (18, "18973553064409449260472376235331"),
            (19, "4536630338860581369328873910626665"),
            (20, "1151178454966303268991128664243557042"),
            (21, "308991125227760514842992561654679405221"),
            (22, "87470525099250663833460093841873159882770"),
            (23, "26045634993717076980553312324382165496411343"),
            (24, "8138039298777944270381420460637129863949889849"),
            (25, "2662347418559335512464065752229073742895672945088"),
        ];

        for (n, val) in expected {
            assert_eq!(
                count_hypercube_nets(n),
                Integer::from_str(val).unwrap(),
                "n={}",
                n
            );
        }
    }
}
