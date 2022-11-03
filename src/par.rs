use crate::scalar::Scalar;
use crate::solve::*;
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

    rhs.par_chunks_exact_mut(n).try_for_each_with(
        vec![S::zero(); n],
        |mut work, mut b| -> Result<(), String> {
            if !trans {
                lsolve(
                    n,
                    &lu.lu_nz,
                    &lu.lu_row_ind,
                    &lu.l_col_ptr,
                    &lu.u_col_ptr,
                    &lu.row_perm,
                    &lu.col_perm,
                    &b,
                    &mut work,
                )?;
                usolve(
                    n,
                    &lu.lu_nz,
                    &lu.lu_row_ind,
                    &lu.l_col_ptr,
                    &lu.u_col_ptr,
                    &lu.row_perm,
                    &lu.col_perm,
                    &mut work,
                    &mut b,
                )?;
            } else {
                utsolve(
                    n,
                    &lu.lu_nz,
                    &lu.lu_row_ind,
                    &lu.l_col_ptr,
                    &lu.u_col_ptr,
                    &lu.row_perm,
                    &lu.col_perm,
                    &b,
                    &mut work,
                )?;
                ltsolve(
                    n,
                    &lu.lu_nz,
                    &lu.lu_row_ind,
                    &lu.l_col_ptr,
                    &lu.u_col_ptr,
                    &lu.row_perm,
                    &lu.col_perm,
                    &mut work,
                    &mut b,
                )?;
            }
            return Ok(());
        },
    )
}
