use crate::internal::{dordstat, OFF};

use crate::rlu::PivotPolicy;

pub fn lucopy(
    pivot: &PivotPolicy,
    pthresh: f64,
    dthresh: f64,
    nzcount: usize,
    jcol1: usize,
    _ncol: usize,
    lastlu: &mut usize,
    lu: &mut [f64],
    lurow: &mut [isize],
    lcolst: &mut [usize],
    ucolst: &mut [usize],
    rperm: &mut [usize],
    cperm: &[usize],
    dense: &mut [f64],
    pattern: &[usize],
    twork: &mut [f64],
) -> Result<isize, String> {
    // Local variables:
    //   nzptr       Index into lurow of current nonzero.
    //   nzst, nzend Loop bounds for nzptr.
    //   irow        Row number of current nonzero (according to A, not PA).
    //   pivrow      Pivot row number (according to A, not PA).
    //   maxpiv      Temporary to find maximum element in column for pivoting.
    //   utemp       Temporary for computing maxpiv.
    //   ujj         Diagonal element U(jcol,jcol) = PtU(pivrow,jcol).
    //   ujjptr      Index into lu and lurow of diagonal element.
    //   dptr        Temporary index into lu and lurow.
    //   diagptr     Index to diagonal element of QAQt
    //   diagpiv     Value of diagonal element

    let jcol = jcol1 - 1; // zero based column

    // Copy column jcol from dense to sparse, recording the position of
    // the diagonal element.
    let mut ujjptr: usize = 0;

    if *pivot == PivotPolicy::NoPivoting || *pivot == PivotPolicy::NoDiagonalElement {
        // No pivoting, diagonal element has irow = jcol.
        // Copy the column elements of U and L, throwing out zeros.

        if ucolst[jcol + 1] - 1 < ucolst[jcol] {
            return Err("zero length (U-I+L) column".to_string());
        }

        // Start with U.
        let mut nzcpy = ucolst[jcol] - 1;
        for nzptr in ucolst[jcol] - 1..lcolst[jcol] - 1 {
            let irow = lurow[nzptr] as usize - 1;

            if pattern[irow] != 0 || irow == cperm[jcol] - 1 {
                lurow[nzcpy] = irow as isize + 1;
                lu[nzcpy] = dense[irow];
                dense[irow] = 0.0;
                nzcpy += 1;
            } else {
                dense[irow] = 0.0;
            }
        }
        let lastu = nzcpy;

        // Now do L. Same action as U, except that we search for diagonal.
        for nzptr in lcolst[jcol] - 1..ucolst[jcol + 1] - 1 {
            let irow = lurow[nzptr] as usize - 1;
            //if irow == cperm[jcol-off] {
            if pattern[irow] == 2 {
                ujjptr = nzcpy + 1
            }
            if pattern[irow] != 0 || /*irow == cperm(jcol-off)) then*/ pattern[irow] == 2 {
                lurow[nzcpy] = irow as isize + 1;
                lu[nzcpy] = dense[irow];
                dense[irow] = 0.0;
                nzcpy += 1;
            } else {
                dense[irow] = 0.0;
            }
        }

        lcolst[jcol] = lastu + 1;
        ucolst[jcol + 1] = nzcpy + 1;
        *lastlu = nzcpy;

        if *pivot == PivotPolicy::NoDiagonalElement {
            let zpivot = 0; //pivrow
            return Ok(zpivot);
        }
    } else {
        let udthreshabs: f64;
        let ldthreshabs: f64;
        let mut rnd: usize = 0;

        // Partial and threshold pivoting.
        if ucolst[jcol + 1] - 1 < lcolst[jcol] {
            return Err("zero length L column".to_string());
        }

        // Partial pivoting, diagonal elt. has max. magnitude in L.
        // Compute the drop threshold for the column
        if nzcount <= 0 {
            let mut maxpivglb = -1.0;
            for nzptr in ucolst[jcol] - 1..lcolst[jcol] - 1 {
                let irow = lurow[nzptr] as usize;
                let utemp = dense[irow - OFF].abs();
                if utemp > maxpivglb {
                    maxpivglb = utemp;
                }
            }
            udthreshabs = dthresh * maxpivglb;

            maxpivglb = -1.0;
            for nzptr in lcolst[jcol] - 1..ucolst[jcol + 1] - 1 {
                let irow = lurow[nzptr] as usize;
                let utemp = dense[irow - OFF].abs();
                if utemp > maxpivglb {
                    maxpivglb = utemp;
                }
            }
            ldthreshabs = dthresh * maxpivglb;
        } else {
            let mut i: usize = 0;
            for nzptr in ucolst[jcol] - 1..lcolst[jcol] - 1 {
                let irow = lurow[nzptr] as usize;
                twork[i] = dense[irow - OFF].abs();
                i += 1;
            }
            if nzcount < i {
                let mut kth: f64 = 0.0;
                let mut info: isize = 0;
                dordstat(&mut rnd, i, i - nzcount + 1, twork, &mut kth, &mut info);
                udthreshabs = kth;
            } else {
                udthreshabs = 0.0;
            }

            let mut i: usize = 0;
            for nzptr in lcolst[jcol] - 1..ucolst[jcol + 1] - 1 {
                let irow = lurow[nzptr] as usize;
                twork[i] = dense[irow - OFF].abs();
                i += 1;
            }
            if nzcount < i {
                let mut kth: f64 = 0.0;
                let mut info: isize = 0;
                dordstat(&mut rnd, i, i - nzcount + 1, twork, &mut kth, &mut info);
                ldthreshabs = kth;
            } else {
                ldthreshabs = 0.0;
            }
        }

        // Copy the column elements of U, throwing out zeros.
        let mut nzcpy = ucolst[jcol] - 1;
        if lcolst[jcol] - 1 >= ucolst[jcol] {
            for nzptr in ucolst[jcol] - 1..lcolst[jcol] - 1 {
                let irow = lurow[nzptr] as usize - 1;

                //if (pattern(irow) .ne. 0 .or. pattern(irow) .eq. 2) then
                if pattern[irow] != 0 || dense[irow].abs() >= udthreshabs {
                    lurow[nzcpy] = irow as isize + 1;
                    lu[nzcpy] = dense[irow];
                    dense[irow] = 0.0;
                    nzcpy += 1;
                } else {
                    dense[irow] = 0.0;
                }
            }
        }
        let lastu = nzcpy;

        // Copy the column elements of L, throwing out zeros.
        // Keep track of maximum magnitude element for pivot.

        if ucolst[jcol + 1] - 1 < lcolst[jcol] {
            return Err("zero length L column".to_string());
        }

        // Partial pivoting, diagonal elt. has max. magnitude in L.
        let mut diagptr: usize = 0;
        let mut diagpiv: f64 = 0.0;

        ujjptr = 0;
        let mut maxpiv = -1.0;
        let mut maxpivglb = -1.0;

        for nzptr in lcolst[jcol] - 1..ucolst[jcol + 1] - 1 {
            let irow = lurow[nzptr] as usize - 1;
            let utemp = dense[irow].abs();

            //if irow == cperm[jcol-off] {
            if pattern[irow] == 2 {
                diagptr = irow + 1;
                diagpiv = utemp;
                //if diagpiv == 0 { print*, 'WARNING: Numerically zero diagonal element at col', jcol }
            }

            // original
            //if utemp > maxpiv {

            // do not pivot outside the pattern
            // if utemp > maxpiv && pattern[irow] != 0 {

            // Pivot outside pattern.
            if utemp > maxpiv {
                ujjptr = irow + 1;
                maxpiv = utemp;
            }

            // Global pivot outside pattern.
            if utemp > maxpivglb {
                maxpivglb = utemp;
            }
        }

        // Threshold pivoting.
        if diagptr != 0 && diagpiv >= (pthresh * maxpiv) {
            ujjptr = diagptr
        }

        if diagptr == 0 && ujjptr == 0 {
            print!("error: {}", ucolst[jcol + 1] - lcolst[jcol]);
        }

        //if diagptr != ujjptr {
        //	print("pivoting", pthresh, maxpiv, diagpiv, diagptr)
        //}

        diagptr = ujjptr;
        ujjptr = 0;

        for nzptr in lcolst[jcol] - 1..ucolst[jcol + 1] - 1 {
            let irow = lurow[nzptr] as usize - 1;
            let utemp = dense[irow].abs();

            //if irow == cperm[jcol-off] {
            //   diagptr = nzcpy
            //   diagpiv = utemp
            //}

            //if utemp > maxpiv {
            //   ujjptr = nzcpy
            //   maxpiv = utemp
            //}

            // Pattern dropping.
            // if pattern[irow] == 0 && irow != diagptr - 1 {

            // Pattern + threshold dropping.

            if pattern[irow] == 0 && irow != diagptr - 1 && utemp < ldthreshabs {
                dense[irow] = 0.0;
            } else {
                if irow == diagptr - 1 {
                    ujjptr = nzcpy + 1;
                }

                lurow[nzcpy] = irow as isize + 1;
                lu[nzcpy] = dense[irow];
                dense[irow] = 0.0;
                nzcpy += 1;
            }
        }

        lcolst[jcol] = lastu + 1;
        ucolst[jcol + 1] = nzcpy + 1;
        *lastlu = nzcpy;
    }

    // Diagonal element has been found. Swap U(jcol,jcol) from L into U.

    if ujjptr == 0 {
        return Err(format!(
            "ujjptr not set (1) {} {} {}", /*diagptr*/
            ujjptr,
            lcolst[jcol],
            ucolst[jcol + 1] - 1
        ));
    }

    let pivrow = lurow[ujjptr - OFF];
    let ujj = lu[ujjptr - OFF];

    if ujj == 0.0 {
        return Err(format!(
            "numerically zero diagonal element at column {}",
            jcol1
        ));
    }
    let dptr = lcolst[jcol];
    lurow[ujjptr - OFF] = lurow[dptr - OFF];
    lu[ujjptr - OFF] = lu[dptr - OFF];
    lurow[dptr - OFF] = pivrow;
    lu[dptr - OFF] = ujj;
    lcolst[jcol] = dptr + 1;

    // Record the pivot in P.

    rperm[pivrow as usize - OFF] = jcol1;
    //if pivrow == 38 {
    //	print("exchanging", jcol, pivrow)
    //}

    // Divide column jcol of L by U(jcol,jcol).

    let nzst = lcolst[jcol];
    let nzend = ucolst[jcol + 1] - 1;
    if nzst > nzend {
        let zpivot = pivrow;
        return Ok(zpivot);
    }
    for nzptr in nzst..=nzend {
        lu[nzptr - OFF] = lu[nzptr - OFF] / ujj;
    }

    let zpivot = pivrow;
    return Ok(zpivot);
}
