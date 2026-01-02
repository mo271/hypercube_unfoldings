use std::cell::RefCell;

use rayon::{iter::repeat, prelude::*};
use rug::{
    Assign, Integer,
    ops::{AddFrom, Pow, SubFrom},
};
use thread_local::ThreadLocal;

struct PrecomputedData {
    factorials: Vec<Integer>,
    pows_of_evens_times_factorials: Vec<Vec<Integer>>,
    t0_partial_a0_1s: Vec<Integer>,
    num_partitions_with_max: Vec<Vec<u64>>,
    num_partitions: Vec<u64>,
}

impl PrecomputedData {
    fn new(n: usize) -> Self {
        let factorials: Vec<Integer> = (0..=n).map(|i| factorial(i)).collect();
        let t0_partial_a0_1s: Vec<Integer> = (0..=n as u32)
            .map(|x| {
                if x >= 2 {
                    Integer::from(2 * x - 2).pow(x) * Integer::from(2 * x).pow(x - 2)
                } else {
                    Integer::ZERO
                }
            })
            .collect();
        let pows_of_evens_times_factorials: Vec<Vec<Integer>> = (0..=n)
            .map(|x| {
                (0..=n)
                    .map(|b| Integer::from(x * 2).pow(b as u32) * &factorials[b])
                    .collect()
            })
            .collect();

        let mut num_partitions_with_max = vec![vec![0u64; n + 1]; n + 1];
        for i in 0..=n {
            num_partitions_with_max[0][i] = 1;
        }
        for i in 1..=n {
            for k in 1..=i {
                num_partitions_with_max[i][k] =
                    num_partitions_with_max[i][k - 1] + num_partitions_with_max[i - k][k];
            }
            for k in i + 1..=n {
                num_partitions_with_max[i][k] = num_partitions_with_max[i][k - 1];
            }
        }

        let num_partitions = (0..=n).map(|i| num_partitions_with_max[i][i]).collect();
        Self {
            factorials,
            pows_of_evens_times_factorials,
            t0_partial_a0_1s,
            num_partitions_with_max,
            num_partitions,
        }
    }

    fn t0_partial_a0_1(&self, p: &[u32]) -> &Integer {
        // [0] is 0.
        &self.t0_partial_a0_1s[p.get(1).copied().unwrap_or_default() as usize]
    }

    fn unrank_partition(&self, n: u32, rank: u64) -> Vec<u32> {
        let n = n as usize;
        debug_assert!(rank < self.num_partitions_with_max[n][n]);
        let mut ans = vec![0; n + 1];
        let mut current_n = n;
        let mut max_val = n;
        let mut index = rank;
        while current_n > 0 {
            for k in (1..=max_val).rev() {
                let count = current_n
                    .checked_sub(k)
                    .map(|x| self.num_partitions_with_max[x][k])
                    .unwrap_or_default();
                if index < count {
                    ans[k] += 1;
                    current_n -= k;
                    max_val = k;
                } else {
                    index -= count;
                }
            }
        }
        ans
    }
}

pub fn count_hypercube_nets(n: u32) -> Integer {
    let precomputed = PrecomputedData::new(n as usize);

    let mut sums = ThreadLocal::<RefCell<Integer>>::new();

    let bn_size = Integer::from(2u32).pow(n) * &precomputed.factorials[n as usize];
    (0..=n)
        .into_par_iter()
        .flat_map(|n_p| {
            let n_m = n - n_p;
            let num_p = precomputed.num_partitions[n_p as usize] as usize;
            let num_m = precomputed.num_partitions[n_m as usize] as usize;
            let num = num_p.checked_mul(num_m).expect("too many terms");
            repeat(n_p).zip(0..num)
        })
        .for_each(|(n_p, i)| {
            let n_m = n - n_p;
            let num_m = precomputed.num_partitions[n_m as usize] as usize;
            let p_idx = i / num_m;
            let m_idx = i % num_m;
            let p = precomputed.unrank_partition(n_p, p_idx as u64);
            let m = precomputed.unrank_partition(n_m, m_idx as u64);
            let sum = sums.get_or_default();
            sum.borrow_mut()
                .add_from(calculate_term(&p, &m, n, &bn_size, &precomputed));
        });
    let mut tot_sum = Integer::ZERO;
    for i in sums.iter_mut() {
        tot_sum.add_from(&*i.borrow());
    }
    tot_sum / bn_size
}

fn calculate_term(
    p: &[u32],
    m: &[u32],
    n: u32,
    bn_size: &Integer,
    precomputed_data: &PrecomputedData,
) -> Integer {
    let p2 = p.get(2).copied().unwrap_or_default();
    let m1 = m.get(1).copied().unwrap_or_default();
    let a_0: u32;
    if p.len() > 1 && p[1] > 0 {
        a_0 = 1;
    } else if p2 > 0 || m1 > 0 {
        a_0 = 2;
    } else {
        return Integer::ZERO;
    }

    let t_0 = if a_0 == 1 {
        precomputed_data.t0_partial_a0_1(p).clone()
    } else {
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

    // w cannot overflow a u64 even with n = 1000
    let mut w = vec![0u64; (2 * n + 1) as usize];

    for b in 1..p.len() {
        if p[b] > 0 {
            for a in ((2 * b)..w.len()).step_by(b) {
                w[a] += 2 * b as u64 * p[b] as u64;
            }
        }
    }
    for k in 1..m.len() {
        if m[k] > 0 {
            let b = 2 * k;
            for a in ((2 * b)..w.len()).step_by(b) {
                w[a] += b as u64 * m[k] as u64;
            }
        }
    }

    let mut prod = Integer::ONE.clone();

    for a in (a_0 + 1)..=(2 * n) {
        let a_idx = a as usize;
        let p_a = p.get(a_idx).copied().unwrap_or_default();
        let a_idx_h = if a_idx % 2 == 0 {
            a_idx / 2
        } else {
            usize::MAX
        };
        let m_a = m.get(a_idx_h).copied().unwrap_or_default();
        let n_a = p_a * 2 + m_a;
        if n_a == 0 {
            continue;
        }

        let w_a = w[a_idx];
        let term_add = a as u64 * n_a as u64;
        let base1 = term_add + w_a;
        let exp1 = n_a - 1 - p_a;

        let base2 = base1.saturating_sub(2);
        let exp2 = p_a;

        let term = w_a * Integer::from(base1).pow(exp1) * Integer::from(base2).pow(exp2);
        prod *= term;
    }

    let tr_g = t_0 * prod;

    let mut cent_size = Integer::ONE.clone();
    for k in 1..p.len() {
        if p[k] > 0 {
            cent_size *= &precomputed_data.pows_of_evens_times_factorials[k][p[k] as usize];
        }
    }
    for k in 1..m.len() {
        if m[k] > 0 {
            cent_size *= &precomputed_data.pows_of_evens_times_factorials[k][m[k] as usize];
        }
    }

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

fn factorial(n: usize) -> Integer {
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
