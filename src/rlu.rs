#[derive(Debug, PartialEq)]
pub enum PivotPolicy {
    NoDiagonalElement,
    NoPivoting,
    PartialPivoting,
    ThresholdPivoting,
}

pub struct Options {
    pub pivotPolicy: PivotPolicy,
    pub pivotThreshold: f64,
    pub dropThreshold: f64,
    pub colFillRatio: f64,
    pub fillRatio: f64,
    pub expandRatio: f64,
    pub colPerm: Option<Vec<usize>>,
}

impl Options {
    pub fn new() -> Options {
        Options {
            pivotPolicy: PivotPolicy::PartialPivoting,
            pivotThreshold: 1.0,
            dropThreshold: 0.0, // do not drop
            colFillRatio: -1.0, // do not limit column fill ratio
            fillRatio: 4.0,
            expandRatio: 1.2,
            colPerm: None,
        }
    }
}
