mod comp;
mod copy;
mod dfs;
mod internal;
mod maxmatch;
pub mod rlu;
mod solve;

use comp::lucomp;
use copy::lucopy;
use dfs::ludfs;
use internal::ifill;
use maxmatch::maxmatch;
use rlu::*;
use solve::*;
use std::cmp::min;

/// LU is a lower-upper numeric factorization.
#[derive(Debug)]
pub struct LU {
    luSize: usize,
    luNZ: Vec<f64>,
    luRowInd: Vec<isize>,
    lColPtr: Vec<usize>,
    uColPtr: Vec<usize>,

    rowPerm: Vec<usize>,
    colPerm: Vec<usize>,

    nA: usize,
}

/// Factor performs sparse LU factorization with partial pivoting.
///
/// Given a matrix A in sparse format by columns, it performs an LU
/// factorization, with partial or threshold pivoting, if desired. The
/// factorization is PA = LU, where L and U are triangular. P, L, and U
/// are returned.  This subroutine uses the Coleman-Gilbert-Peierls
/// algorithm, in which total time is O(nonzero multiplications).
pub fn factor(
    nA: usize,
    rowind: &[usize],
    colptr: &[usize],
    nzA: &[f64],
    opts: &Options,
) -> Result<LU, String> {
    let nrow = nA;
    let ncol = nA;
    let nnzA = nzA.len();

    if nnzA > nA * nA {
        return Err(format!("nnz ({}) must be < n*n ({})", nnzA, nA * nA));
    }
    if rowind.len() != nzA.len() {
        return Err(format!(
            "len rowind ({}) must be nnz ({})",
            rowind.len(),
            nzA.len()
        ));
    }
    if colptr.len() != ncol + 1 {
        return Err(format!(
            "len colptr ({}) must be ncol+1 ({})",
            colptr.len(),
            ncol + 1
        ));
    }

    match &opts.colPerm {
        Some(colPerm) => {
            // If a column permutation is specified, it must be a length ncol permutation.
            if colPerm.len() != ncol {
                return Err(format!(
                    "column permutation ({}) must be a length ncol {}",
                    colPerm.len(),
                    ncol
                ));
            }
            for v in colPerm.iter() {
                if v < &0 || v >= &ncol {
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
    let mut colptrA = vec![0; nA + 1];
    let mut rowindA = vec![0; nnzA];
    for jcol in 0..nA + 1 {
        colptrA[jcol] = colptr[jcol] + 1;
    }
    for jcol in 0..nnzA {
        rowindA[jcol] = rowind[jcol] + 1;
    }

    // Allocate work arrays.
    let mut rwork = vec![0.0; nrow];
    let mut twork = vec![0.0; nrow];
    let mut found = vec![0; nrow];
    let mut child = vec![0; nrow];
    let mut parent = vec![0; nrow];
    let mut pattern = vec![0; nrow];

    // Create lu structure.
    let luSize = ((nnzA as f64) * opts.fillRatio) as usize;
    let mut lu = LU {
        luSize,
        luNZ: vec![0.0; luSize],
        luRowInd: vec![0; luSize],
        uColPtr: vec![0; ncol + 1],
        lColPtr: vec![0; ncol],
        rowPerm: vec![0; nrow],
        colPerm: vec![0; ncol],
        nA,
    };

    let (mut rmatch, mut cmatch) = maxmatch(
        nrow,
        ncol,
        &colptrA,
        &rowindA,
        &mut lu.lColPtr,
        &mut lu.uColPtr,
        &mut lu.rowPerm,
        &mut lu.colPerm,
        &mut lu.luRowInd,
    )?;

    for jcol in 0..ncol {
        if cmatch[jcol] == 0 {
            println!("warning: perfect matching not found");
            break;
        }
    }

    // Initialize useful values and zero out the dense vectors.
    // If we are threshold pivoting, get row counts.
    let mut lastlu = 0;

    let mut localPivotPolicy = &opts.pivotPolicy;
    //let lasta = colptrA[ncol] - 1;
    lu.uColPtr[0] = 1;

    //ifill(pattern, nrow, 0)
    //ifill(found, nrow, 0)
    //rfill(rwork, nrow, 0)
    ifill(&mut lu.rowPerm, nrow, 0);

    match &opts.colPerm {
        Some(colPerm) => {
            for jcol in 0..ncol {
                lu.colPerm[jcol] = colPerm[jcol] + 1;
            }
        }
        None => {
            for jcol in 0..ncol {
                lu.colPerm[jcol] = jcol + 1
            }
        }
    }

    // Compute one column at a time.
    for jcol in 1..=ncol {
        // Mark pointer to new column, ensure it is large enough.
        if lastlu + nrow >= lu.luSize {
            let newSize = ((lu.luSize as f64) * opts.expandRatio) as usize;

            println!("expanding LU to {} nonzeros", newSize);

            let mut luNZ = vec![0.0; newSize];
            luNZ[..lu.luSize].copy_from_slice(&lu.luNZ[..]);
            lu.luNZ = luNZ;

            let mut luRowInd = vec![0; newSize];
            luRowInd[..lu.luSize].copy_from_slice(&lu.luRowInd[..]);
            lu.luRowInd = luRowInd;

            lu.luSize = newSize;
        }

        // Set up nonzero pattern.
        let (origRow, thisCol) = {
            let jjj = lu.colPerm[jcol - 1];
            for i in colptrA[jjj - 1]..colptrA[jjj] {
                pattern[rowindA[i - 1] - 1] = 1;
            }

            let thisCol = lu.colPerm[jcol - 1];
            let origRow = cmatch[thisCol - 1];

            pattern[origRow - 1] = 2;

            if lu.rowPerm[origRow - 1] != 0 {
                return Err("pivot row from max-matching already used".to_string());
            }

            (origRow, thisCol)
        };

        // Depth-first search from each above-diagonal nonzero of column
        // jcol of A, allocating storage for column jcol of U in
        // topological order and also for the non-fill part of column
        // jcol of L.
        ludfs(
            jcol,
            nzA,
            &rowindA,
            &colptrA,
            &mut lastlu,
            &mut lu.luRowInd,
            &mut lu.lColPtr,
            &lu.uColPtr,
            &lu.rowPerm,
            &lu.colPerm,
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
            &lu.luNZ,
            &mut lu.luRowInd,
            &lu.lColPtr,
            &mut lu.uColPtr,
            &lu.rowPerm,
            &lu.colPerm,
            &mut rwork,
            &mut found,
            &pattern,
        );

        // Copy the dense vector into the sparse data structure, find the
        // diagonal element (pivoting if specified), and divide the
        // column of L by it.
        let nzCountLimit =
            (opts.colFillRatio * ((colptrA[thisCol] - colptrA[thisCol - 1] + 1) as f64)) as usize;

        let zpivot = lucopy(
            &localPivotPolicy,
            opts.pivotThreshold,
            opts.dropThreshold,
            nzCountLimit,
            jcol,
            ncol,
            &mut lastlu,
            &mut lu.luNZ,
            &mut lu.luRowInd,
            &mut lu.lColPtr,
            &mut lu.uColPtr,
            &mut lu.rowPerm,
            &mut lu.colPerm,
            &mut rwork,
            &mut pattern,
            &mut twork,
        )?;
        if zpivot == -1 {
            return Err(format!("lucopy: jcol={}", jcol));
        }

        {
            let jjj = lu.colPerm[jcol - 1];
            for i in colptrA[jjj - 1]..colptrA[jjj] {
                pattern[rowindA[i - 1] - 1] = 0;
            }

            pattern[origRow - 1] = 0;

            let pivtRow = zpivot;
            let othrCol = rmatch[pivtRow as usize - 1];

            cmatch[thisCol - 1] = pivtRow as usize;
            cmatch[othrCol - 1] = origRow;
            rmatch[origRow - 1] = othrCol;
            rmatch[pivtRow as usize - 1] = thisCol;

            //pattern[thisCol - 1] = 0
        }

        // If there are no diagonal elements after this column, change the pivot mode.
        if jcol == nrow {
            localPivotPolicy = &PivotPolicy::NoDiagonalElement;
        }
    }

    // Fill in the zero entries of the permutation vector, and renumber the
    // rows so the data structure represents L and U, not PtL and PtU.
    let mut jcol = ncol + 1;
    for i in 0..nrow {
        if lu.rowPerm[i] == 0 {
            lu.rowPerm[i] = jcol;
            jcol = jcol + 1;
        }
    }

    for i in 0..lastlu {
        lu.luRowInd[i] = lu.rowPerm[lu.luRowInd[i] as usize - 1] as isize;
    }

    Ok(lu)
}

/// Solve Ax=b for one or more right-hand-sides given the numeric
/// factorization of A from `factor`.
pub fn solve(lu: &LU, rhs: &mut [&mut [f64]], trans: bool) -> Result<(), String> {
    let n = lu.nA;
    if rhs.len() == 0 {
        return Err("one or more rhs must be specified".to_string());
    }
    for (i, b) in rhs.iter().enumerate() {
        if b.len() != n {
            return Err(format!(
                "len b[{}] ({}) must equal ord(A) ({})",
                i,
                b.len(),
                n
            ));
        }
    }
    let mut work = vec![0.0; n];

    for i in 0..rhs.len() {
        let mut b = &mut rhs[i];
        if !trans {
            lsolve(
                n,
                &lu.luNZ,
                &lu.luRowInd,
                &lu.lColPtr,
                &lu.uColPtr,
                &lu.rowPerm,
                &lu.colPerm,
                &b,
                &mut work,
            )?;
            usolve(
                n,
                &lu.luNZ,
                &lu.luRowInd,
                &lu.lColPtr,
                &lu.uColPtr,
                &lu.rowPerm,
                &lu.colPerm,
                &mut work,
                &mut b,
            )?;
        } else {
            utsolve(
                n,
                &lu.luNZ,
                &lu.luRowInd,
                &lu.lColPtr,
                &lu.uColPtr,
                &lu.rowPerm,
                &lu.colPerm,
                &b,
                &mut work,
            )?;
            ltsolve(
                n,
                &lu.luNZ,
                &lu.luRowInd,
                &lu.lColPtr,
                &lu.uColPtr,
                &lu.rowPerm,
                &lu.colPerm,
                &mut work,
                &mut b,
            )?;
        }
    }

    return Ok(());
}
