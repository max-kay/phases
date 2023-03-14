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

impl RandAtom for BinAtoms {
    type Concentration = f64;

    fn uniform(rng: &mut Pcg64) -> Self {
        if rng.gen() {
            Self::A
        } else {
            Self::B
        }
    }

    fn with_concentration(rng: &mut Pcg64, c_a: Self::Concentration) -> Self {
        if rng.gen_bool(c_a) {
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
