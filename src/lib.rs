use std::cell::RefCell;

use indicatif::{ProgressBar, ProgressStyle};
use rayon::{iter::repeat, prelude::*};
use rug::{
    Assign, Integer,
    ops::{AddFrom, MulFrom, Pow, SubFrom},
};
use thread_local::ThreadLocal;

struct PrecomputedData {
    factorials: Vec<Integer>,
    pows_of_evens_times_factorials: Vec<Vec<Integer>>,
    t0_partial_a0_1s: Vec<Integer>,
    num_partitions_with_max: Vec<Vec<u64>>,
    num_partitions: Vec<u64>,
}

#[derive(Default)]
struct ThreadCache {
    matrices: Vec<Vec<Vec<Integer>>>,
    w: Vec<u64>,
    t_0: Integer,
    cent: Integer,
    val: Integer,
    val_tmp: Integer,
    count: u64,
}

impl PrecomputedData {
    fn new(n: usize) -> Self {
        let factorials: Vec<Integer> = (0..=n).map(factorial).collect();
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
        num_partitions_with_max[0].fill(1);
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

    fn unrank_partition(&self, n: usize, rank: u64, ans: &mut Vec<u32>) {
        debug_assert!(rank < self.num_partitions_with_max[n][n]);
        ans.resize(n + 1, 0);
        ans.fill(0);
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
    }
}

pub fn count_hypercube_nets(n: u32) -> Integer {
    let n = n as usize;
    let precomputed = PrecomputedData::new(n);

    let mut sums = ThreadLocal::<RefCell<Integer>>::new();
    let mut caches = ThreadLocal::<RefCell<ThreadCache>>::new();
    let unrank_cache_p = ThreadLocal::<RefCell<Vec<u32>>>::new();
    let unrank_cache_m = ThreadLocal::<RefCell<Vec<u32>>>::new();

    let total_iterations: u64 = (0..=n)
        .map(|n_p| {
            let n_m = n - n_p;
            precomputed.num_partitions[n_p] * precomputed.num_partitions[n_m]
        })
        .sum();

    let pb = if total_iterations > 1000 {
        let pb = ProgressBar::new(total_iterations);
        pb.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
                )
                .unwrap()
                .progress_chars("#>-"),
        );
        pb
    } else {
        ProgressBar::hidden()
    };

    let bn_size = Integer::from(2u32).pow(n as u32) * &precomputed.factorials[n];
    (0..=n)
        .into_par_iter()
        .flat_map(|n_p| {
            let n_m = n - n_p;
            let num_p = precomputed.num_partitions[n_p] as usize;
            let num_m = precomputed.num_partitions[n_m] as usize;
            let num = num_p.checked_mul(num_m).expect("too many terms");
            repeat(n_p).zip(0..num)
        })
        .for_each(|(n_p, i)| {
            let n_m = n - n_p;
            let num_m = precomputed.num_partitions[n_m] as usize;
            let p_idx = i / num_m;
            let m_idx = i % num_m;
            let p = unrank_cache_p.get_or_default();
            let m = unrank_cache_m.get_or_default();
            let mut p = p.borrow_mut();
            let mut m = m.borrow_mut();
            precomputed.unrank_partition(n_p, p_idx as u64, &mut p);
            precomputed.unrank_partition(n_m, m_idx as u64, &mut m);
            let sum = sums.get_or_default();
            let cache = caches.get_or_default();
            let mut cache = cache.borrow_mut();
            calculate_term(
                &p,
                &m,
                n,
                &bn_size,
                &precomputed,
                &mut *sum.borrow_mut(),
                &mut *cache,
            );
            cache.count += 1;
            if cache.count >= 1000 {
                pb.inc(1000);
                cache.count -= 1000;
            }
        });

    let mut tot_sum = Integer::ZERO;
    for i in sums.iter_mut() {
        tot_sum.add_from(&*i.borrow());
    }
    for cache in caches.iter_mut() {
        pb.inc(cache.borrow().count);
    }
    pb.finish_and_clear();
    tot_sum / bn_size
}

fn calculate_term(
    p: &[u32],
    m: &[u32],
    n: usize,
    bn_size: &Integer,
    precomputed_data: &PrecomputedData,
    accumulator: &mut Integer,
    cache: &mut ThreadCache,
) {
    let ThreadCache {
        matrices,
        w,
        t_0,
        cent,
        val,
        val_tmp,
        ..
    } = cache;
    let p2 = p.get(2).copied().unwrap_or_default();
    let m1 = m.get(1).copied().unwrap_or_default();
    let a_0: usize;
    if p.len() > 1 && p[1] > 0 {
        a_0 = 1;
    } else if p2 > 0 || m1 > 0 {
        a_0 = 2;
    } else {
        return;
    }

    if a_0 == 1 {
        t_0.assign(precomputed_data.t0_partial_a0_1(p))
    } else {
        let v = (2 * p2 + m1) as usize;

        if v == 0 || p2 == 0 {
            return;
        }

        let size = v - 1;
        if matrices.len() <= size {
            matrices.resize(size + 1, vec![]);
        }
        let mat = &mut matrices[size];
        if mat.len() != size {
            *mat = vec![vec![Integer::ZERO; size]; size];
        }
        let num_paired_rows = 2 * p2 as usize;

        for (r, row) in mat.iter_mut().enumerate() {
            for (c, val) in row.iter_mut().enumerate() {
                let orig_r = r + 1;
                let orig_c = c + 1;
                if r == c {
                    if orig_r < num_paired_rows {
                        val.assign(2 * v as i64 - 3);
                    } else {
                        val.assign(2 * v as i64 - 2);
                    }
                } else {
                    let is_pair = if orig_r < num_paired_rows && orig_c < num_paired_rows {
                        (orig_r ^ 1) == orig_c
                    } else {
                        false
                    };
                    if is_pair {
                        val.assign(-1);
                    } else {
                        val.assign(-2);
                    }
                }
            }
        }
        let det = determinant_bigint(mat, val, val_tmp);
        if !det.is_positive() {
            return;
        }
        t_0.assign(2 * p2);
        t_0.mul_from(det);
    };

    if t_0.is_zero() {
        return;
    }

    // w cannot overflow a u64 even with n = 1000
    if w.len() != 2 * n + 1 {
        w.resize(2 * n + 1, 0);
    }
    w.fill(0);

    for (b, &v) in p.iter().enumerate().skip(1) {
        if v > 0 {
            for a in ((2 * b)..w.len()).step_by(b) {
                w[a] += 2 * b as u64 * v as u64;
            }
        }
    }
    for (k, &v) in m.iter().enumerate().skip(1) {
        if v > 0 {
            let b = 2 * k;
            for a in ((2 * b)..w.len()).step_by(b) {
                w[a] += b as u64 * v as u64;
            }
        }
    }

    for a in (a_0 + 1)..=(2 * n) {
        let a_idx = a;
        let p_a = p.get(a_idx).copied().unwrap_or_default();
        let a_idx_h = if a_idx.is_multiple_of(2) {
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
        *t_0 *= term;
    }

    cent.assign(1);
    for (k, &v) in p.iter().enumerate().skip(1) {
        if v > 0 {
            *cent *= &precomputed_data.pows_of_evens_times_factorials[k][v as usize];
        }
    }
    for (k, &v) in m.iter().enumerate().skip(1) {
        if v > 0 {
            *cent *= &precomputed_data.pows_of_evens_times_factorials[k][v as usize];
        }
    }

    *t_0 *= bn_size;
    t_0.div_exact_mut(cent);
    accumulator.add_from(&*t_0);
}

fn determinant_bigint<'a>(
    matrix: &'a mut [Vec<Integer>],
    val: &mut Integer,
    val_tmp: &mut Integer,
) -> &'a Integer {
    let n = matrix.len();
    if n == 0 {
        return Integer::ONE;
    }
    if n == 1 {
        return &matrix[0][0];
    }

    for k in 0..n - 1 {
        if matrix[k][k].is_zero() {
            let mut swap_row = None;
            for (i, row) in matrix.iter_mut().enumerate().skip(k) {
                if !row[k].is_zero() {
                    swap_row = Some(i);
                    break;
                }
            }
            if let Some(r) = swap_row {
                matrix.swap(k, r);
                matrix[k].iter_mut().for_each(|v| {
                    v.sub_from(&Integer::ZERO);
                });
            } else {
                matrix[0][0].assign(0);
                return &matrix[0][0];
            }
        }
        let (before, after) = matrix.split_at_mut(k);
        let denom = if k == 0 {
            Integer::ONE
        } else {
            &before[k - 1][k - 1]
        };
        if denom.is_zero() {
            matrix[0][0].assign(0);
            return &matrix[0][0];
        }
        for i in k + 1..n {
            for j in k + 1..n {
                val_tmp.assign(&after[i - k][j] * &after[0][k]);
                val.assign(&after[i - k][k] * &after[0][j]);
                // val = val_tmp - val
                val.sub_from(&*val_tmp);
                after[i - k][j].assign(val.div_exact_ref(denom));
            }
        }
    }
    &matrix[n - 1][n - 1]
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
