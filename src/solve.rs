use crate::internal::OFF;

pub fn lsolve(
    n: usize,
    lu: &[f64],
    lurow: &[isize],
    lcolst: &[usize],
    ucolst: &[usize],
    rperm: &[usize],
    _cperm: &[usize],
    b: &[f64],
    x: &mut [f64],
) -> Result<(), String> {
    if n <= 0 {
        return Err(format!("lsolve called with nonpositive n = {}", n));
    }
    // Solve the system.
    for i in 1..=n {
        x[rperm[i - OFF] - OFF] = b[i - OFF];
    }
    for j in 1..=n {
        let nzst = lcolst[j - OFF];
        let nzend = ucolst[j + 1 - OFF] - 1;
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
            let i = lurow[nzptr - OFF] as usize;
            if i <= j || i > n {
                return Err(format!(
                    "lsolve, illegal row i in column j of L: i={}, j={}, nzptr={}",
                    i, j, nzptr
                ));
            }
            x[i - OFF] -= lu[nzptr - OFF] * x[j - OFF];
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
    _cperm: &[usize],
    b: &mut [f64],
    x: &mut [f64],
) -> Result<(), String> {
    if n <= 0 {
        return Err(format!("ltsolve called with nonpositive n={}", n));
    }
    // Solve the system.
    for i in 1..=n {
        x[i - OFF] = b[i - OFF];
    }
    for j in (1..=n).rev() {
        let nzst = lcolst[j - OFF];
        let nzend = ucolst[j + 1 - OFF] - 1;
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
            let i = lurow[nzptr - OFF] as usize;
            if i <= j || i > n {
                return Err(format!(
                    "ltsolve, illegal row i in column j of L: i={}, j={}, nzptr={}",
                    i, j, nzptr
                ));
            }
            x[j - OFF] -= lu[nzptr - OFF] * x[i - OFF];
        }
    }

    for i in 1..=n {
        b[i - OFF] = x[i - OFF];
    }

    for i in 1..=n {
        //x[rperm[i-off]-off] = b[i-off]
        x[i - OFF] = b[rperm[i - OFF] - OFF];
    }

    Ok(())
}

pub fn usolve(
    n: usize,
    lu: &[f64],
    lurow: &[isize],
    lcolst: &[usize],
    ucolst: &[usize],
    _rperm: &[usize],
    cperm: &[usize],
    b: &mut [f64],
    x: &mut [f64],
) -> Result<(), String> {
    if n <= 0 {
        return Err(format!("usolve called with nonpositive n={}", n));
    }
    for i in 1..=n {
        x[i - OFF] = b[i - OFF];
    }

    for jj in 1..=n {
        let j = n + 1 - jj;
        let nzst = ucolst[j - OFF];
        let mut nzend = lcolst[j - OFF] - 1;
        if nzst < 1 || nzst > nzend {
            return Err(format!(
                "usolve, inconsistent column of U: j={}, nzst={}, nzend={}",
                j, nzst, nzend
            ));
        }
        if lurow[nzend - OFF] != j as isize {
            return Err(format!("usolve, diagonal elt of col j is not in last place: j={}, nzend={}, lurow[nzend]={}", j, nzend, lurow[nzend- OFF]));
        }
        if lu[nzend - OFF] == 0.0 {
            return Err(format!("usolve, zero diagonal element in column j={}", j));
        }
        x[j - OFF] = x[j - OFF] / lu[nzend - OFF];
        nzend = nzend - 1;
        if nzst > nzend {
            continue;
        }
        for nzptr in nzst..=nzend {
            let i = lurow[nzptr - OFF] as usize;
            if i <= 0 || i >= j {
                return Err(format!(
                    "usolve, illegal row i in column j of U: i={}, j={}, nzptr={}",
                    i, j, nzptr
                ));
            }
            x[i - OFF] -= lu[nzptr - OFF] * x[j - OFF];
        }
    }

    for i in 1..=n {
        b[i - OFF] = x[i - OFF];
    }
    for i in 1..=n {
        x[cperm[i - OFF] - OFF] = b[i - OFF];
    }

    Ok(())
}

pub fn utsolve(
    n: usize,
    lu: &[f64],
    lurow: &[isize],
    lcolst: &[usize],
    ucolst: &[usize],
    _rperm: &[usize],
    cperm: &[usize],
    b: &[f64],
    x: &mut [f64],
) -> Result<(), String> {
    if n <= 0 {
        return Err(format!("utsolve called with nonpositive n={}", n));
    }
    for i in 1..=n {
        x[i - OFF] = b[cperm[i - OFF] - OFF];
    }

    for j in 1..=n {
        let nzst = ucolst[j - OFF];
        let mut nzend = lcolst[j - OFF] - 1;
        if nzst < 1 || nzst > nzend {
            return Err(format!(
                "utsolve, inconsistent column of U: j={}, nzst={}, nzend={}",
                j, nzst, nzend
            ));
        }
        if lurow[nzend - OFF] != j as isize {
            return Err(format!("utsolve, diagonal elt of col j is not in last place: j={}, nzend={}, lurow[nzend]={}", j, nzend, lurow[nzend- OFF]));
        }
        if lu[nzend - OFF] == 0.0 {
            return Err(format!("utsolve, zero diagonal element in column j={}", j));
        }
        nzend = nzend - 1;
        if nzst > nzend {
            x[j - OFF] = x[j - OFF] / lu[nzend + 1 - OFF];
            continue;
        }
        for nzptr in nzst..=nzend {
            let i = lurow[nzptr - OFF] as usize;
            if i <= 0 || i >= j {
                return Err(format!(
                    "utsolve, illegal row i in column j of U: i={}, j={}, nzptr={}",
                    i, j, nzptr
                ));
            }
            x[j - OFF] -= lu[nzptr - OFF] * x[i - OFF];
        }
        x[j - OFF] = x[j - OFF] / lu[nzend + 1 - OFF];
    }

    Ok(())
}
