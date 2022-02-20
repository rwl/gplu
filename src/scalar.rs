use num_complex::Complex;
use num_traits::Zero;
use std::ops::{Div, Mul, SubAssign};

pub trait Scalar<Output = Self>:
    PartialEq + Clone + Mul<Output = Output> + Div<Output = Output> + SubAssign + Copy + Zero
{
    fn abs(self) -> f64;
}

impl Scalar for f64 {
    fn abs(self) -> f64 {
        self.abs()
    }
}

impl Scalar for f32 {
    fn abs(self) -> f64 {
        self.abs() as f64
    }
}

impl Scalar for Complex<f64> {
    fn abs(self) -> f64 {
        self.norm()
    }
}

impl Scalar for Complex<f32> {
    fn abs(self) -> f64 {
        self.norm() as f64
    }
}
