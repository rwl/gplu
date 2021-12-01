pub const off: usize = 1;

// ifill fills an integer array with a given value.
pub fn ifill(a: &mut [usize], la: usize, ival: usize) {
    for i in 1..=la {
        a[i - off] = ival;
    }
}

pub fn dordstat(
    rnd: &mut usize,
    n: usize,
    k: usize,
    A: &mut [f64],
    kth: &mut f64,
    info: &mut isize,
) {
    let mut i: usize;
    let mut j: usize;
    let mut x: f64;

    if k < 0 || k > n {
        *info = -1;
        return;
    }

    let mut p = 1;
    let mut r = n;

    loop {
        if p == r {
            break;
        }

        if r - p >= 8 {
            *rnd = (1366 * *rnd + 150889) % 714025;
            let q = p + (*rnd % (r - p + 1));

            let tmp = A[p - off];
            A[p - off] = A[q - off];
            A[q - off] = tmp;
        }

        x = A[p - off];
        i = p - 1;
        j = r + 1;

        loop {
            j = j - 1;
            while A[j - off] > x {
                j = j - 1;
            }

            i = i + 1;
            while A[i - off] < x {
                i = i + 1;
            }

            if i < j {
                let tmp = A[i - off];
                A[i - off] = A[j - off];
                A[j - off] = tmp;
                continue;
            }

            if j < k {
                p = j + 1;
            } else {
                r = j;
            }
            break;
        }
    }

    *kth = A[p - off];
    *info = 0;
}
