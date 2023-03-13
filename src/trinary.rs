use rand::Rng;
use rand_pcg::Pcg64;

use crate::{ModularGrid, RandAtom};

impl From<&ModularGrid<TriAtoms>> for &[u8] {
    fn from(value: &ModularGrid<TriAtoms>) -> Self {
        unsafe {
            std::slice::from_raw_parts(
                value.grid.as_ptr() as *const u8,
                value.grid.len() * std::mem::size_of::<TriAtoms>(),
            )
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
#[repr(u8)]
pub enum TriAtoms {
    #[default]
    A,
    B,
    C,
}

impl RandAtom for TriAtoms {
    type Concentration = (f64, f64);

    fn uniform(rng: &mut Pcg64) -> Self {
        match rng.gen_range(0..3){
            0 => Self::A,
            1 => Self::B,
            2 => Self::C,
            _ => unreachable!()
        }
    }

    fn with_concentration(rng: &mut Pcg64, c_a: Self::Concentration) -> Self {
        if rng.gen_bool(c_a.0) {
            Self::A
        } else if rng.gen_bool(c_a.1 / (1.0 - c_a.0)) {
            Self::B
        } else {
            Self::C
        }
    }
}
