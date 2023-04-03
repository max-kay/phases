use rand::Rng;
use rand_pcg::Pcg64;

use crate::{ModularArray, RandAtom};

#[derive(Debug, Clone, Copy, Default)]
#[repr(u8)]
pub enum TerAtoms {
    #[default]
    A,
    B,
    C,
}

#[derive(Debug, Clone, Copy)]
pub struct TerConcentration {
    c_a: f64,
    c_b_without_a: f64,
}

impl TerConcentration {
    pub fn new(x1: f64, x2: f64, x3: f64) -> Self {
        Self {
            c_a: x1 / (x1 + x2 + x3),
            c_b_without_a: x2 / (x2 + x3),
        }
    }
}

impl RandAtom for TerAtoms {
    type Concentration = TerConcentration;

    fn uniform(rng: &mut Pcg64) -> Self {
        match rng.gen_range(0..3) {
            0 => Self::A,
            1 => Self::B,
            2 => Self::C,
            _ => unreachable!(),
        }
    }

    fn with_concentration(rng: &mut Pcg64, c: Self::Concentration) -> Self {
        if rng.gen_bool(c.c_a) {
            Self::A
        } else if rng.gen_bool(c.c_b_without_a) {
            Self::B
        } else {
            Self::C
        }
    }
}

impl<const WIDTH: usize, const HEIGHT: usize> From<&ModularArray<TerAtoms, WIDTH, HEIGHT>>
    for &[u8]
where
    [(); WIDTH * HEIGHT]:,
{
    fn from(value: &ModularArray<TerAtoms, WIDTH, HEIGHT>) -> Self {
        unsafe {
            std::slice::from_raw_parts(
                value.grid.as_ptr() as *const u8,
                value.grid.len() * std::mem::size_of::<TerAtoms>(),
            )
        }
    }
}
