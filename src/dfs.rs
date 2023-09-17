use crate::internal::OFF;
use crate::scalar::Scalar;

/// Depth-first search to allocate storage for U.
///
/// # Input parameters
///
/// ```txt
///   jcol             current column number.
///   a, arow, acolst  the matrix A; see lufact for format.
///   rperm            row permutation P.
///                    perm(r) = s > 0 means row r of A is row s < jcol of PA.
///                    perm(r) = 0 means row r of A has not yet been used as a
///                    pivot and is therefore still below the diagonal.
///   cperm            column permutation.
/// ```
///
/// # Modified parameters (see below for exit values):
///
/// ```txt
///   lastlu           last used position in lurow array.
///   lurow, lcolst, ucolst  nonzero structure of Pt(L-I+U);
///                          see LU for format.
///   dense            current column as a dense vector.
///   found            integer array for marking nonzeros in this column of
///                    Pt(L-I+U) that have been allocated space in lurow.
///                    Also, marks reached columns in depth-first search.
///                    found(i)=jcol if i was found in this column.
///   parent           parent(i) is the parent of vertex i in the dfs,
///                    or 0 if i is a root of the search.
///   child            child(i) is the index in lurow of the next unexplored
///                    child of vertex i.
///                    Note that parent and child are also indexed according to
///                    the vertex numbering of A, not PA; thus child(i) is
///                    the position of a nonzero in column rperm(i),
///                    not column i.
/// ```
///
/// # On entry
///
/// ```txt
///   found(*)<jcol
///   dense(*)=0.0
///   ucolst(jcol)=lastlu+1 is the first free index in lurow.
/// ```
///
/// # On exit
///
/// ```txt
///   found(i)=jcol iff i is a nonzero of column jcol of PtU or
///     a non-fill nonzero of column jcol of PtL.
///   dense(*)=column jcol of A.
///     Note that found and dense are kept according to the row
///     numbering of A, not PA.
///   lurow has the rows of the above-diagonal nonzeros of col jcol of U in
///     reverse topological order, followed by the non-fill nonzeros of col
///     jcol of L and the diagonal elt of U, in no particular order.
///     These rows also are numbered according to A, not PA.
///   lcolst(jcol) is the index of the first nonzero in col j of L.
///   lastlu is the index of the last non-fill nonzero in col j of L.
/// ```
pub fn ludfs<S: Scalar>(
    jcol: usize,
    a: &[S],
    arow: &[usize],
    acolst: &[usize],
    lastlu: &mut usize,
    lurow: &mut [isize],
    lcolst: &mut [usize],
    ucolst: &[usize],
    rperm: &[usize],
    cperm: &[usize],
    dense: &mut [S],
    found: &mut [usize],
    parent: &mut [usize],
    child: &mut [usize],
) -> Result<(), String> {
    // Depth-first search through columns of L from each nonzero of
    // column jcol of A that is above the diagonal in PA.

    // For each krow such that A(krow,jcol) is nonzero do...

    // Range of indices in arow for column jcol of A.
    let nzast = acolst[cperm[jcol - OFF] - OFF];
    let nzaend = acolst[cperm[jcol - OFF]]; //+1-off]

    if nzaend < nzast {
        return Err(format!(
            "ludfs, negative length for column {} of A. nzast={} nzend={}",
            jcol, nzast, nzaend
        ));
    }
    let nzaend = nzaend - 1;
    for nzaptr in nzast - 1..nzaend {
        // pointer to current position in arow (zero based)
        // Current vertex in depth-first search (numbered according to A, not PA) (zero based).
        let mut krow: isize = arow[nzaptr] as isize - 1;

        // Copy A(krow,jcol) into the dense vector. If above diagonal in
        // PA, start a depth-first search in column rperm(krow) of L.

        dense[krow as usize] = a[nzaptr];
        if rperm[krow as usize] == 0
            || found[krow as usize] == jcol
            || dense[krow as usize] == S::zero()
        {
            continue;
        }
        parent[krow as usize] = 0;
        found[krow as usize] = jcol;
        let mut chdptr = lcolst[rperm[krow as usize] - OFF]; // Index of current child of current vertex.

        // The main depth-first search loop starts here.
        // repeat
        //   if krow has a child that is not yet found
        //   then step forward
        //   else step back
        // until a step back leads to 0
        'l100: loop {
            // Look for an unfound child of krow.
            let chdend = ucolst[rperm[krow as usize]]; // Next index after last child of current vertex.

            'l200: loop {
                if chdptr < chdend {
                    // Possible next vertex in depth-first search (zero based).
                    let nextk: isize = lurow[chdptr - OFF] - 1;
                    chdptr = chdptr + 1;
                    if rperm[nextk as usize] == 0 {
                        continue 'l200;
                    }
                    if found[nextk as usize] == jcol {
                        continue 'l200;
                    }
                    // Take a step forward.

                    //l300:
                    child[krow as usize] = chdptr;
                    parent[nextk as usize] = krow as usize + 1;
                    krow = nextk;
                    found[krow as usize] = jcol;
                    chdptr = lcolst[rperm[krow as usize] - OFF];
                    //goto l100
                    continue 'l100;
                }
                break;
            }
            // Take a step back.

            // Allocate space for U(rperm(k),jcol) = PtU(krow,jcol) in the sparse data structure.
            *lastlu = *lastlu + 1;
            lurow[*lastlu - OFF] = krow + 1;
            krow = parent[krow as usize] as isize - 1;
            if krow >= 0 {
                chdptr = child[krow as usize];
                //goto l100
                continue 'l100;
            }
            // The main depth-first search loop ends here.
            break;
        }
    }
    // Close off column jcol of U and allocate space for the non-fill
    // entries of column jcol of L.
    // The diagonal element goes in L, not U, until we do the column
    // division at the end of the major step.

    lcolst[jcol - OFF] = *lastlu + 1;
    for nzaptr in nzast - 1..nzaend {
        let krow = arow[nzaptr];
        if rperm[krow - OFF] == 0 {
            found[krow - OFF] = jcol;
            *lastlu += 1;
            lurow[*lastlu - OFF] = krow as isize;
        }
    }

    Ok(())
}
