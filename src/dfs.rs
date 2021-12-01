use crate::internal::off;

pub fn ludfs(
    jcol: usize,
    a: &[f64],
    arow: &[usize],
    acolst: &[usize],
    lastlu: &mut usize,
    lurow: &mut [isize],
    lcolst: &mut [usize],
    ucolst: &[usize],
    rperm: &[usize],
    cperm: &[usize],
    dense: &mut [f64],
    found: &mut [usize],
    parent: &mut [usize],
    child: &mut [usize],
) -> Result<(), String> {
    // Depth-first search through columns of L from each nonzero of
    // column jcol of A that is above the diagonal in PA.

    // For each krow such that A(krow,jcol) is nonzero do...

    // Range of indices in arow for column jcol of A.
    let nzast = acolst[cperm[jcol - off] - off];
    let nzaend = acolst[cperm[jcol - off]]; //+1-off]

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
        if rperm[krow as usize] == 0 || found[krow as usize] == jcol || dense[krow as usize] == 0.0
        {
            continue;
        }
        parent[krow as usize] = 0;
        found[krow as usize] = jcol;
        let mut chdptr = lcolst[rperm[krow as usize] - off]; // Index of current child of current vertex.

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
                    let nextk: isize = lurow[chdptr - off] - 1;
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
                    chdptr = lcolst[rperm[krow as usize] - off];
                    //goto l100
                    continue 'l100;
                }
                break;
            }
            // Take a step back.

            // Allocate space for U(rperm(k),jcol) = PtU(krow,jcol) in the sparse data structure.
            *lastlu = *lastlu + 1;
            lurow[*lastlu - off] = krow + 1;
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

    lcolst[jcol - off] = *lastlu + 1;
    for nzaptr in nzast - 1..nzaend {
        let krow = arow[nzaptr];
        if rperm[krow - off] == 0 {
            found[krow - off] = jcol;
            *lastlu += 1;
            lurow[*lastlu - off] = krow as isize;
        }
    }

    Ok(())
}
