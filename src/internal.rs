pub const OFF: usize = 1;

// ifill fills an integer array with a given value.
pub fn ifill(a: &mut [usize], la: usize, ival: usize) {
    for i in 1..=la {
        a[i - OFF] = ival;
    }
}

pub fn dordstat(
    rnd: &mut usize,
    n: usize,
    k: usize,
    a: &mut [f64],
    kth: &mut f64,
    info: &mut isize,
) {
    let mut i: usize;
    let mut j: usize;
    let mut x: f64;

    if
    /*k < 0 ||*/
    k > n {
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

            let tmp = a[p - OFF];
            a[p - OFF] = a[q - OFF];
            a[q - OFF] = tmp;
        }

        x = a[p - OFF];
        i = p - 1;
        j = r + 1;

        loop {
            j = j - 1;
            while a[j - OFF] > x {
                j = j - 1;
            }

            i = i + 1;
            while a[i - OFF] < x {
                i = i + 1;
            }

            if i < j {
                let tmp = a[i - OFF];
                a[i - OFF] = a[j - OFF];
                a[j - OFF] = tmp;
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

    *kth = a[p - OFF];
    *info = 0;
}

#[cfg(feature = "debug")]
macro_rules! debug_println {
    ($( $args:expr ),*) => { println!( $( $args ),* ); }
}

#[cfg(not(feature = "debug"))]
macro_rules! debug_println {
    ($( $args:expr ),*) => {};
}

pub(crate) use debug_println;
