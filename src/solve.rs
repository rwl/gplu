use crate::scalar::Scalar;
use crate::LU;

pub fn lu_solve<S: Scalar>(
    lu: &LU<S>,
    b: &mut [S],
    work: &mut [S],
    trans: bool,
) -> Result<(), String> {
    if b.len() != lu.n {
        return Err(format!("len b ({}) must be n ({})", b.len(), lu.n));
    }

    if !trans {
        lsolve(
            lu.n,
            &lu.lu_nz,
            &lu.lu_row_ind,
            &lu.l_col_ptr,
            &lu.u_col_ptr,
            &lu.row_perm,
            &lu.col_perm,
            b,
            work,
        )?;
        usolve(
            lu.n,
            &lu.lu_nz,
            &lu.lu_row_ind,
            &lu.l_col_ptr,
            &lu.u_col_ptr,
            &lu.row_perm,
            &lu.col_perm,
            work,
            b,
        )?;
    } else {
        utsolve(
            lu.n,
            &lu.lu_nz,
            &lu.lu_row_ind,
            &lu.l_col_ptr,
            &lu.u_col_ptr,
            &lu.row_perm,
            &lu.col_perm,
            b,
            work,
        )?;
        ltsolve(
            lu.n,
            &lu.lu_nz,
            &lu.lu_row_ind,
            &lu.l_col_ptr,
            &lu.u_col_ptr,
            &lu.row_perm,
            &lu.col_perm,
            work,
            b,
        )?;
    }

    return Ok(());
}

pub fn lsolve<S: Scalar>(
    n: usize,
    lu: &[S],
    lurow: &[isize],
    lcolst: &[usize],
    ucolst: &[usize],
    rperm: &[usize],
    _cperm: &[usize],
    b: &[S],
    x: &mut [S],
) -> Result<(), String> {
    #[cfg(feature = "debug")]
    if n == 0 {
        return Err(format!("lsolve called with nonpositive n = {}", n));
    }
    // Solve the system.
    for i in 0..n {
        x[rperm[i] - 1] = b[i];
    }

    for j in 0..n {
        let nzst = lcolst[j] - 1;
        let nzend = ucolst[j + 1] - 1;
        #[cfg(feature = "debug")]
        if nzst > nzend {
            return Err(format!(
                "lsolve, inconsistent column of L: j={} nzst={}, nzend={}",
                j, nzst, nzend
            ));
        }
        for nzptr in nzst..nzend {
            let i = (lurow[nzptr] - 1) as usize;
            #[cfg(feature = "debug")]
            if i <= j || i >= n {
                return Err(format!(
                    "lsolve, illegal row i in column j of L: i={}, j={}, nzptr={}",
                    i, j, nzptr
                ));
            }
            x[i] -= lu[nzptr] * x[j];
        }
    }
    Ok(())
}

pub fn ltsolve<S: Scalar>(
    n: usize,
    lu: &[S],
    lurow: &[isize],
    lcolst: &[usize],
    ucolst: &[usize],
    rperm: &[usize],
    _cperm: &[usize],
    b: &mut [S],
    x: &mut [S],
) -> Result<(), String> {
    #[cfg(feature = "debug")]
    if n == 0 {
        return Err(format!("ltsolve called with nonpositive n={}", n));
    }
    // Solve the system.
    x.copy_from_slice(b);

    for j in (0..n).rev() {
        let nzst = lcolst[j] - 1;
        let nzend = ucolst[j + 1] - 1;
        #[cfg(feature = "debug")]
        if nzst > nzend {
            return Err(format!(
                "ltsolve, inconsistent column of L: j={}, nzst={}, nzend={}",
                j, nzst, nzend
            ));
        }
        for nzptr in nzst..nzend {
            let i = (lurow[nzptr] as usize) - 1;
            #[cfg(feature = "debug")]
            if i <= j || i >= n {
                return Err(format!(
                    "ltsolve, illegal row i in column j of L: i={}, j={}, nzptr={}",
                    i, j, nzptr
                ));
            }
            x[j] -= lu[nzptr] * x[i];
        }
    }

    b.copy_from_slice(x);

    for i in 0..n {
        x[i] = b[rperm[i] - 1];
    }

    Ok(())
}

pub fn usolve<S: Scalar>(
    n: usize,
    lu: &[S],
    lurow: &[isize],
    lcolst: &[usize],
    ucolst: &[usize],
    _rperm: &[usize],
    cperm: &[usize],
    b: &mut [S],
    x: &mut [S],
) -> Result<(), String> {
    #[cfg(feature = "debug")]
    if n == 0 {
        return Err(format!("usolve called with nonpositive n={}", n));
    }
    x.copy_from_slice(b);

    for j in (0..n).rev() {
        let nzst = ucolst[j] - 1;
        let nzend = lcolst[j] - 1;
        #[cfg(feature = "debug")]
        if nzst >= nzend {
            return Err(format!(
                "usolve, inconsistent column of U: j={}, nzst={}, nzend={}",
                j, nzst, nzend
            ));
        }
        #[cfg(feature = "debug")]
        if lurow[nzend - 1] - 1 != j as isize {
            return Err(format!("usolve, diagonal elt of col j is not in last place: j={}, nzend={}, lurow[nzend]={}", j, nzend, lurow[nzend - 1]));
        }
        #[cfg(feature = "debug")]
        if lu[nzend - 1] == S::zero() {
            return Err(format!("usolve, zero diagonal element in column j={}", j));
        }
        x[j] = x[j] / lu[nzend - 1];

        let nzend = nzend - 1;

        for nzptr in nzst..nzend {
            let i = (lurow[nzptr] - 1) as usize;
            #[cfg(feature = "debug")]
            if lurow[nzptr] <= 0 || i >= j {
                return Err(format!(
                    "usolve, illegal row i in column j of U: i={}, j={}, nzptr={}",
                    i, j, nzptr
                ));
            }
            x[i] -= lu[nzptr] * x[j];
        }
    }

    b.copy_from_slice(x);

    for i in 0..n {
        x[cperm[i] - 1] = b[i];
    }

    Ok(())
}

pub fn utsolve<S: Scalar>(
    n: usize,
    lu: &[S],
    lurow: &[isize],
    lcolst: &[usize],
    ucolst: &[usize],
    _rperm: &[usize],
    cperm: &[usize],
    b: &[S],
    x: &mut [S],
) -> Result<(), String> {
    #[cfg(feature = "debug")]
    if n == 0 {
        return Err(format!("utsolve called with nonpositive n={}", n));
    }
    for i in 0..n {
        x[i] = b[cperm[i] - 1];
    }

    for j in 0..n {
        let nzst = ucolst[j] - 1;
        let nzend = lcolst[j] - 1;
        #[cfg(feature = "debug")]
        if nzst >= nzend {
            return Err(format!(
                "utsolve, inconsistent column of U: j={}, nzst={}, nzend={}",
                j, nzst, nzend
            ));
        }
        #[cfg(feature = "debug")]
        if lurow[nzend - 1] - 1 != j as isize {
            return Err(format!("utsolve, diagonal elt of col j is not in last place: j={}, nzend={}, lurow[nzend]={}", j, nzend, lurow[nzend - 1]));
        }
        #[cfg(feature = "debug")]
        if lu[nzend - 1] == S::zero() {
            return Err(format!("utsolve, zero diagonal element in column j={}", j));
        }
        x[j] = x[j] / lu[nzend - 1];

        let nzend = nzend - 1;

        for nzptr in nzst..nzend {
            let i = (lurow[nzptr] - 1) as usize;
            #[cfg(feature = "debug")]
            if lurow[nzptr] <= 0 || i >= j {
                return Err(format!(
                    "utsolve, illegal row i in column j of U: i={}, j={}, nzptr={}",
                    i, j, nzptr
                ));
            }
            x[j] -= lu[nzptr] * x[i];
        }
    }

    Ok(())
}
