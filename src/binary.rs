use rand::Rng;
use rand_pcg::Pcg64;

use crate::{ModularArray, RandAtom};

#[derive(Debug, Clone, Copy, Default)]
#[repr(u8)]
pub enum BinAtoms {
    #[default]
    A,
    B,
}

#[derive(Debug, Clone, Copy)]
pub struct BinConcentration {
    c_a: f64,
}

impl BinConcentration {
    pub fn new(x1: f64, x2: f64) -> Self {
        Self {
            c_a: x1 / (x1 + x2),
        }
    }
}

impl RandAtom for BinAtoms {
    type Concentration = BinConcentration;

    fn uniform(rng: &mut Pcg64) -> Self {
        if rng.gen() {
            Self::A
        } else {
            Self::B
        }
    }

    fn with_concentration(rng: &mut Pcg64, c: Self::Concentration) -> Self {
        if rng.gen_bool(c.c_a) {
            Self::A
        } else {
            Self::B
        }
    }
}

impl<const WIDTH: usize, const HEIGHT: usize> From<&ModularArray<BinAtoms, WIDTH, HEIGHT>>
    for &[u8]
where
    [(); WIDTH * HEIGHT]:,
{
    fn from(value: &ModularArray<BinAtoms, WIDTH, HEIGHT>) -> Self {
        unsafe {
            std::slice::from_raw_parts(
                value.grid.as_ptr() as *const u8,
                value.grid.len() * std::mem::size_of::<BinAtoms>(),
            )
        }
    }
}
