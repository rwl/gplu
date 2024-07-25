#[derive(Debug, PartialEq)]
pub enum PivotPolicy {
    NoDiagonalElement,
    NoPivoting,
    PartialPivoting,
    ThresholdPivoting,
}

pub struct Options {
    pub pivot_policy: PivotPolicy,
    /// The fraction of max pivot candidate acceptable for pivoting.
    /// Default value is 1.
    pub pivot_threshold: f64,

    /// For each major step of the algorithm, the pivot is chosen to
    /// be a nonzero below the diagonal in the current column of A
    /// with the most nonzeros to the right in its row, with absolute
    /// value at least drop_threshold*maxpiv, where maxpiv is the
    /// largest absolute value below the diagonal in the current column.
    /// Note that if drop_threshold <= 0.0, then the pivot is chosen
    /// purely on the basis of row sparsity. Also, if
    /// drop_threshold >= 1.0, then the pivoting is effectively partial
    /// pivoting with ties broken on the basis of sparsity.
    pub drop_threshold: f64,

    /// Column fill ratio. If < 0 the column fill ratio is not limited.
    /// Default value is -1.
    pub col_fill_ratio: f64,

    /// sets the ratio of the initial LU size to NNZ.
    /// Default value is 4.
    pub fill_ratio: f64,

    /// Sets the ratio for LU size growth.
    pub expand_ratio: f64,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            pivot_policy: PivotPolicy::PartialPivoting,
            pivot_threshold: 1.0,
            drop_threshold: 0.0,  // do not drop
            col_fill_ratio: -1.0, // do not limit column fill ratio
            fill_ratio: 4.0,
            expand_ratio: 1.2,
        }
    }
}
