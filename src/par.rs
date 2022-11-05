use crate::scalar::Scalar;
use crate::solve::lu_solve;
use crate::LU;
use rayon::iter::ParallelIterator;
use rayon::slice::ParallelSliceMut;

/// Solve `Ax=b` for one or more right-hand-sides given the numeric
/// factorization of A from `factor`.
pub fn par_solve<S: Scalar + Send + Sync>(
    lu: &LU<S>,
    rhs: &mut [S],
    trans: bool,
) -> Result<(), String> {
    let n = lu.n;
    if rhs.len() % n != 0 {
        return Err(format!(
            "len rhs ({}) must be a multiple of n ({})",
            rhs.len(),
            n
        ));
    }

    rhs.par_chunks_exact_mut(n)
        .try_for_each_with(vec![S::zero(); n], |work, b| -> Result<(), String> {
            lu_solve(lu, b, work, trans)
        })
}
