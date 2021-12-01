use crate::internal::off;

pub fn lsolve(
    n: usize,
    lu: &[f64],
    lurow: &[isize],
    lcolst: &[usize],
    ucolst: &[usize],
    rperm: &[usize],
    cperm: &[usize],
    b: &[f64],
    x: &mut [f64],
) -> Result<(), String> {
    if n <= 0 {
        return Err(format!("lsolve called with nonpositive n = {}", n));
    }
    // Solve the system.
    for i in 1..=n {
        x[rperm[i - off] - off] = b[i - off];
    }
    for j in 1..=n {
        let nzst = lcolst[j - off];
        let nzend = ucolst[j + 1 - off] - 1;
        if nzst < 1 || nzst > nzend + 1 {
            return Err(format!(
                "lsolve, inconsistent column of L: j={} nzst={}, nzend={}",
                j, nzst, nzend
            ));
        }
        if nzst > nzend {
            continue;
        }
        for nzptr in nzst..=nzend {
            let i = lurow[nzptr - off] as usize;
            if i <= j || i > n {
                return Err(format!(
                    "lsolve, illegal row i in column j of L: i={}, j={}, nzptr={}",
                    i, j, nzptr
                ));
            }
            x[i - off] -= lu[nzptr - off] * x[j - off];
        }
    }
    Ok(())
}

pub fn ltsolve(
    n: usize,
    lu: &[f64],
    lurow: &[isize],
    lcolst: &[usize],
    ucolst: &[usize],
    rperm: &[usize],
    cperm: &[usize],
    b: &mut [f64],
    x: &mut [f64],
) -> Result<(), String> {
    if n <= 0 {
        return Err(format!("ltsolve called with nonpositive n={}", n));
    }
    // Solve the system.
    for i in 1..=n {
        x[i - off] = b[i - off];
    }
    for j in (n..=1).rev() {
        let nzst = lcolst[j - off];
        let nzend = ucolst[j + 1 - off] - 1;
        if nzst < 1 || nzst > nzend + 1 {
            return Err(format!(
                "ltsolve, inconsistent column of L: j={}, nzst={}, nzend={}",
                j, nzst, nzend
            ));
        }
        if nzst > nzend {
            continue;
        }
        for nzptr in nzst..=nzend {
            let i = lurow[nzptr - off] as usize;
            if i <= j || i > n {
                return Err(format!(
                    "ltsolve, illegal row i in column j of L: i={}, j={}, nzptr={}",
                    i, j, nzptr
                ));
            }
            x[j - off] -= lu[nzptr - off] * x[i - off];
        }
    }

    for i in 1..=n {
        b[i - off] = x[i - off];
    }

    for i in 1..=n {
        //x[rperm[i-off]-off] = b[i-off]
        x[i - off] = b[rperm[i - off] - off];
    }

    Ok(())
}

pub fn usolve(
    n: usize,
    lu: &[f64],
    lurow: &[isize],
    lcolst: &[usize],
    ucolst: &[usize],
    rperm: &[usize],
    cperm: &[usize],
    b: &mut [f64],
    x: &mut [f64],
) -> Result<(), String> {
    if n <= 0 {
        return Err(format!("usolve called with nonpositive n={}", n));
    }
    for i in 1..=n {
        x[i - off] = b[i - off];
    }

    for jj in 1..=n {
        let j = n + 1 - jj;
        let nzst = ucolst[j - off];
        let mut nzend = lcolst[j - off] - 1;
        if nzst < 1 || nzst > nzend {
            return Err(format!(
                "usolve, inconsistent column of U: j={}, nzst={}, nzend={}",
                j, nzst, nzend
            ));
        }
        if lurow[nzend - off] != j as isize {
            return Err(format!("usolve, diagonal elt of col j is not in last place: j={}, nzend={}, lurow[nzend]={}", j, nzend, lurow[nzend-off]));
        }
        if lu[nzend - off] == 0.0 {
            return Err(format!("usolve, zero diagonal element in column j={}", j));
        }
        x[j - off] = x[j - off] / lu[nzend - off];
        nzend = nzend - 1;
        if nzst > nzend {
            continue;
        }
        for nzptr in nzst..=nzend {
            let i = lurow[nzptr - off] as usize;
            if i <= 0 || i >= j {
                return Err(format!(
                    "usolve, illegal row i in column j of U: i={}, j={}, nzptr={}",
                    i, j, nzptr
                ));
            }
            x[i - off] -= lu[nzptr - off] * x[j - off];
        }
    }

    for i in 1..=n {
        b[i - off] = x[i - off];
    }
    for i in 1..=n {
        x[cperm[i - off] - off] = b[i - off];
    }

    Ok(())
}

pub fn utsolve(
    n: usize,
    lu: &[f64],
    lurow: &[isize],
    lcolst: &[usize],
    ucolst: &[usize],
    rperm: &[usize],
    cperm: &[usize],
    b: &[f64],
    x: &mut [f64],
) -> Result<(), String> {
    if n <= 0 {
        return Err(format!("utsolve called with nonpositive n={}", n));
    }
    for i in 1..=n {
        x[i - off] = b[cperm[i - off] - off];
    }

    for j in 1..=n {
        let nzst = ucolst[j - off];
        let mut nzend = lcolst[j - off] - 1;
        if nzst < 1 || nzst > nzend {
            return Err(format!(
                "utsolve, inconsistent column of U: j={}, nzst={}, nzend={}",
                j, nzst, nzend
            ));
        }
        if lurow[nzend - off] != j as isize {
            return Err(format!("utsolve, diagonal elt of col j is not in last place: j={}, nzend={}, lurow[nzend]={}", j, nzend, lurow[nzend-off]));
        }
        if lu[nzend - off] == 0.0 {
            return Err(format!("utsolve, zero diagonal element in column j={}", j));
        }
        nzend = nzend - 1;
        if nzst > nzend {
            x[j - off] = x[j - off] / lu[nzend + 1 - off];
            continue;
        }
        for nzptr in nzst..=nzend {
            let i = lurow[nzptr - off] as usize;
            if i <= 0 || i >= j {
                return Err(format!(
                    "utsolve, illegal row i in column j of U: i={}, j={}, nzptr={}",
                    i, j, nzptr
                ));
            }
            x[j - off] -= lu[nzptr - off] * x[i - off];
        }
        x[j - off] = x[j - off] / lu[nzend + 1 - off];
    }

    Ok(())
}
