//! Implements an algorithm for sparse Gaussian elimination in time
//! proportional to arithmetic operations.
//! It does NOT do any preordering for sparsity. We recommend preordering
//! the columns of the input matrix by using minimum degree on the
//! symmetric structure of `(A-transpose)A`, i.e. the column intersection
//! graph of `A`.
//!
//! ```bibtex
//! @article{Gilbert1988,
//!   doi = {10.1137/0909058},
//!   url = {https://doi.org/10.1137/0909058},
//!   year  = {1988},
//!   month = {sep},
//!   publisher = {Society for Industrial {\&} Applied Mathematics ({SIAM})},
//!   volume = {9},
//!   number = {5},
//!   pages = {862--874},
//!   author = {John R. Gilbert and Tim Peierls},
//!   title = {Sparse Partial Pivoting in Time Proportional to Arithmetic Operations},
//!   journal = {{SIAM} Journal on Scientific and Statistical Computing}
//! }
//! ```

extern crate num_complex;
extern crate num_traits;

mod comp;
mod copy;
mod dfs;
mod internal;
mod maxmatch;
#[cfg(feature = "rayon")]
mod par;
mod rlu;
mod scalar;
mod solve;

pub use crate::rlu::*;
use comp::lucomp;
use copy::lucopy;
use dfs::ludfs;
use internal::*;
use maxmatch::maxmatch;
use num_traits::PrimInt;
#[cfg(feature = "rayon")]
pub use par::par_solve;
pub use scalar::Scalar;
use solve::lu_solve;
use std::fmt::Display;

/// LU is a lower-upper numeric factorization.
#[derive(Debug, Clone)]
pub struct LU<S: Scalar> {
    // Index of last nonzero in lu; number of nonzeros in L-I+U.
    lu_size: usize,
    // Nonzeros in L and U.  Nonzeros in each column are contiguous.
    // Columns of U and L alternate: u1, l1, u2, l2, ..., un, ln.
    // Nonzeros are not necessarily in order by row within columns.
    // The diagonal elements of L, which are all 1, are not stored.
    lu_nz: Vec<S>,
    // lurow(i) is the row number of the nonzero lu(i).
    // During the computation, these correspond to row numbers in A,
    // so we really store the non-triangular PtL and PtU.
    // At the end we transform them to row numbers in PA, so the
    // L and U we return are really triangular.
    lu_row_ind: Vec<isize>,
    // lcolst(j) is the index in lu of the first nonzero in col j of L.
    l_col_ptr: Vec<usize>,
    // ucolst(j) is the index in lu of the first nonzero in col j of U.
    // The last nonzero of col j of U is in position lcolst(j)-1, and
    // the last nonzero of col j of L is in position ucolst(j+1)-1.
    // ucolst(ncol+1) is one more than the number of nonzeros in L-I+U.
    //
    // Notice that ucolst has dimension ncol+1 and lcolst has dimension ncol,
    // although ucolst(ncol+1)=lcolst(ncol) because the last column of L
    // contains only the diagonal one, which is not stored.
    u_col_ptr: Vec<usize>,

    // perm(r) = s means that row r of A is in position s in PA.
    row_perm: Vec<usize>,
    col_perm: Vec<usize>,

    n: usize,
}

/// Factor performs sparse LU factorization with partial pivoting.
///
/// Given a matrix A in sparse format by columns, it performs an LU
/// factorization, with partial or threshold pivoting, if desired. The
/// factorization is `PA = LU`, where `L` and `U` are triangular. `P`, `L`, and `U`
/// are returned.  This subroutine uses the Coleman-Gilbert-Peierls
/// algorithm, in which total time is O(nonzero multiplications).
pub fn factor<I: PrimInt + Display, S: Scalar>(
    nn: I,
    rowind0: &[I],
    colptr0: &[I],
    nz: &[S],
    col_perm: Option<&[I]>,
    opts: &Options,
) -> Result<LU<S>, String> {
    let n = nn.to_usize().unwrap();

    let nrow = n;
    let ncol = n;
    let nnz = nz.len();

    if nnz > n * n {
        return Err(format!("nnz ({}) must be < n*n ({})", nnz, n * n));
    }
    if rowind0.len() != nz.len() {
        return Err(format!(
            "len rowind ({}) must be nnz ({})",
            rowind0.len(),
            nz.len()
        ));
    }
    if colptr0.len() != ncol + 1 {
        return Err(format!(
            "len colptr ({}) must be ncol+1 ({})",
            colptr0.len(),
            ncol + 1
        ));
    }

    match col_perm {
        Some(col_perm) => {
            // If a column permutation is specified, it must be a length ncol permutation.
            if col_perm.len() != ncol {
                return Err(format!(
                    "column permutation ({}) must be a length ncol {}",
                    col_perm.len(),
                    ncol
                ));
            }
            for v in col_perm.iter() {
                if v < &I::zero() || v >= &I::from(ncol).unwrap() {
                    return Err(format!(
                        "column permutation {} out of range [0,{})",
                        v, ncol
                    ));
                }
            }
        }
        None => (),
    }

    // Convert the descriptor to 1-base if necessary.
    let mut colptr = vec![0; n + 1];
    let mut rowind = vec![0; nnz];
    for jcol in 0..n + 1 {
        colptr[jcol] = colptr0[jcol].to_usize().unwrap() + 1;
    }
    for jcol in 0..nnz {
        rowind[jcol] = rowind0[jcol].to_usize().unwrap() + 1;
    }

    // Allocate work arrays.
    let mut rwork = vec![S::zero(); nrow];
    let mut twork = vec![0.0; nrow];
    // integer vector used to control depth-first search.
    let mut found = vec![0; nrow];
    let mut child = vec![0; nrow]; // also used by ludfs.
    let mut parent = vec![0; nrow]; // also used by ludfs.
    let mut pattern = vec![0; nrow];

    // Create lu structure.
    let lu_size = ((nnz as f64) * opts.fill_ratio) as usize;
    let mut lu = LU {
        lu_size: lu_size,
        lu_nz: vec![S::zero(); lu_size],
        lu_row_ind: vec![0; lu_size],
        u_col_ptr: vec![0; ncol + 1],
        l_col_ptr: vec![0; ncol],
        row_perm: vec![0; nrow],
        col_perm: vec![0; ncol],
        n: n,
    };

    let (mut rmatch, mut cmatch) = maxmatch(
        nrow,
        ncol,
        &colptr,
        &rowind,
        &mut lu.l_col_ptr,
        &mut lu.u_col_ptr,
        &mut lu.row_perm,
        &mut lu.col_perm,
        &mut lu.lu_row_ind,
    )?;

    #[cfg(feature = "debug")]
    for jcol in 0..ncol {
        if cmatch[jcol] == 0 {
            debug_println!("warning: perfect matching not found");
            break;
        }
    }

    // Initialize useful values and zero out the dense vectors.
    // If we are threshold pivoting, get row counts.
    let mut lastlu = 0;

    let mut local_pivot_policy = &opts.pivot_policy;
    //let lasta = colptrA[ncol] - 1;
    lu.u_col_ptr[0] = 1;

    //ifill(pattern, nrow, 0)
    //ifill(found, nrow, 0)
    //rfill(rwork, nrow, 0)
    ifill(&mut lu.row_perm, nrow, 0);

    match col_perm {
        Some(col_perm) => {
            for jcol in 0..ncol {
                lu.col_perm[jcol] = col_perm[jcol].to_usize().unwrap() + 1;
            }
        }
        None => {
            for jcol in 0..ncol {
                lu.col_perm[jcol] = jcol + 1
            }
        }
    }

    // Compute one column at a time.
    for jcol in 1..=ncol {
        // Mark pointer to new column, ensure it is large enough.
        if lastlu + nrow >= lu.lu_size {
            let new_size = ((lu.lu_size as f64) * opts.expand_ratio) as usize;

            debug_println!("expanding LU to {} nonzeros", new_size);

            let mut lu_nz = vec![S::zero(); new_size];
            lu_nz[..lu.lu_size].copy_from_slice(&lu.lu_nz[..]);
            lu.lu_nz = lu_nz;

            let mut lu_row_ind = vec![0; new_size];
            lu_row_ind[..lu.lu_size].copy_from_slice(&lu.lu_row_ind[..]);
            lu.lu_row_ind = lu_row_ind;

            lu.lu_size = new_size;
        }

        // Set up nonzero pattern.
        let (orig_row, this_col) = {
            let jjj = lu.col_perm[jcol - 1];
            for i in colptr[jjj - 1]..colptr[jjj] {
                pattern[rowind[i - 1] - 1] = 1;
            }

            let this_col = lu.col_perm[jcol - 1];
            let orig_row = cmatch[this_col - 1];

            pattern[orig_row - 1] = 2;

            if lu.row_perm[orig_row - 1] != 0 {
                return Err("pivot row from max-matching already used".to_string());
            }

            (orig_row, this_col)
        };

        // Depth-first search from each above-diagonal nonzero of column
        // jcol of A, allocating storage for column jcol of U in
        // topological order and also for the non-fill part of column
        // jcol of L.
        ludfs(
            jcol,
            nz,
            &rowind,
            &colptr,
            &mut lastlu,
            &mut lu.lu_row_ind,
            &mut lu.l_col_ptr,
            &lu.u_col_ptr,
            &lu.row_perm,
            &lu.col_perm,
            &mut rwork,
            &mut found,
            &mut parent,
            &mut child,
        )?;

        // Compute the values of column jcol of L and U in the dense
        // vector, allocating storage for fill in L as necessary.

        lucomp(
            jcol,
            &mut lastlu,
            &lu.lu_nz,
            &mut lu.lu_row_ind,
            &lu.l_col_ptr,
            &mut lu.u_col_ptr,
            &lu.row_perm,
            &lu.col_perm,
            &mut rwork,
            &mut found,
            &pattern,
        );

        // Copy the dense vector into the sparse data structure, find the
        // diagonal element (pivoting if specified), and divide the
        // column of L by it.
        let nz_count_limit =
            (opts.col_fill_ratio * ((colptr[this_col] - colptr[this_col - 1] + 1) as f64)) as isize;

        let zpivot = lucopy(
            &local_pivot_policy,
            opts.pivot_threshold,
            opts.drop_threshold,
            nz_count_limit,
            jcol,
            ncol,
            &mut lastlu,
            &mut lu.lu_nz,
            &mut lu.lu_row_ind,
            &mut lu.l_col_ptr,
            &mut lu.u_col_ptr,
            &mut lu.row_perm,
            &mut lu.col_perm,
            &mut rwork,
            &mut pattern,
            &mut twork,
        )?;
        if zpivot == -1 {
            return Err(format!("lucopy: jcol={}", jcol));
        }

        {
            let jjj = lu.col_perm[jcol - 1];
            for i in colptr[jjj - 1]..colptr[jjj] {
                pattern[rowind[i - 1] - 1] = 0;
            }

            pattern[orig_row - 1] = 0;

            let pivt_row = zpivot;
            let othr_col = rmatch[pivt_row as usize - 1];

            cmatch[this_col - 1] = pivt_row as usize;
            cmatch[othr_col - 1] = orig_row;
            rmatch[orig_row - 1] = othr_col;
            rmatch[pivt_row as usize - 1] = this_col;

            //pattern[thisCol - 1] = 0
        }

        // If there are no diagonal elements after this column, change the pivot mode.
        if jcol == nrow {
            local_pivot_policy = &PivotPolicy::NoDiagonalElement;
        }
    }

    // Fill in the zero entries of the permutation vector, and renumber the
    // rows so the data structure represents L and U, not PtL and PtU.
    let mut jcol = ncol + 1;
    for i in 0..nrow {
        if lu.row_perm[i] == 0 {
            lu.row_perm[i] = jcol;
            jcol = jcol + 1;
        }
    }

    for i in 0..lastlu {
        lu.lu_row_ind[i] = lu.row_perm[lu.lu_row_ind[i] as usize - 1] as isize;
    }

    #[cfg(feature = "debug")]
    {
        let mut ujj: f64 = 0.0;
        let mut minujj = f64::INFINITY;

        for jcol in 1..=ncol {
            ujj = (lu.lu_nz[lu.l_col_ptr[jcol - 1] - 2]).abs();
            if ujj < minujj {
                minujj = ujj;
            }
        }

        debug_println!("last = {}, min = {}", ujj, minujj);
    }

    Ok(lu)
}

/// Solve `Ax=b` for one or more right-hand-sides given the numeric
/// factorization of A from `factor`.
pub fn solve<S: Scalar>(lu: &LU<S>, rhs: &mut [S], trans: bool) -> Result<(), String> {
    let n = lu.n;
    if rhs.len() % n != 0 {
        return Err(format!(
            "len rhs ({}) must be a multiple of n ({})",
            rhs.len(),
            n
        ));
    }
    let mut work = vec![S::zero(); n];

    rhs.chunks_exact_mut(n)
        .try_for_each(|b| lu_solve(lu, b, &mut work, trans))
}
