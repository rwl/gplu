use crate::internal::off;

pub fn lucomp(
    jcol: usize,
    lastlu: &mut usize,
    lu: &[f64],
    lurow: &mut [isize],
    lcolst: &[usize],
    ucolst: &mut [usize],
    rperm: &[usize],
    cperm: &[usize],
    dense: &mut [f64],
    found: &mut [usize],
    pattern: &[usize],
) {
    // Local variables:
    //   nzuptr                pointer to current nonzero PtU(krow,jcol).
    //   nzuend, nnzu, nzuind  used to compute nzuptr.
    //   krow                  row index of current nonzero (according to A,
    //                         not PA), which is PtU(krow,jcol) = U(kcol,jcol).
    //   kcol                  rperm(krow), that is, row index of current
    //                         nonzero according to PA.
    //   ukj                   value of PtU(krow,jcol).
    //   nzlptr                pointer to PtL(irow,kcol), being used for update.
    //   nzlst, nzlend         used to compute nzlptr.
    //   irow                  row index in which update is taking place,
    //                         according to A (not PA).

    //    For each krow with PtU(krow,jcol) != 0, in reverse postorder, use
    //    column kcol = rperm(krow) of L to update the current column.
    let nzuend = lcolst[jcol - off];
    let nnzu = nzuend - ucolst[jcol - off];
    if nnzu != 0 {
        for nzuind in 1..=nnzu {
            let nzuptr = nzuend - nzuind;
            let krow = lurow[nzuptr - off] - 1;
            let kcol = rperm[krow as usize] - 1;
            let ukj = dense[krow as usize];
            //if pattern[rperm[krow]-off] == 0 {
            //	ukj = 0
            //}

            // For each irow with PtL(irow,kcol) != 0, update PtL(irow,jcol) or PtU(irow,jcol)

            let nzlst = lcolst[kcol];
            let nzlend = ucolst[kcol + 1] - 1;
            if nzlend < nzlst {
                continue;
            }
            for nzlptr in nzlst - 1..nzlend {
                let irow = lurow[nzlptr] - 1;
                dense[irow as usize] -= ukj * lu[nzlptr];

                // If this is a new nonzero in L, allocate storage for it.
                if found[irow as usize] != jcol {
                    found[irow as usize] = jcol;
                    lurow[*lastlu] = irow + 1;
                    *lastlu += 1
                }
            }
        }
    }
    ucolst[jcol + 1 - off] = *lastlu + 1;
}
