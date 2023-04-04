use rand::Rng;
use rand_pcg::Pcg64;

use crate::{Concentration, ModularArray, RandAtom};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[repr(u8)]
pub enum BinAtoms {
    #[default]
    A,
    B,
}

#[derive(Debug, Clone, Copy)]
pub struct BinConcentration {
    c_a: f64,
    c_b: f64,
}

impl BinConcentration {
    pub fn new(n_a: f64, n_b: f64) -> Self {
        let n_tot = n_a + n_b;
        Self {
            c_a: n_a / n_tot,
            c_b: n_a / n_tot,
        }
    }
}

impl Concentration for BinConcentration {
    fn uniform() -> Self {
        Self::new(1.0, 1.0)
    }

    fn max_entropy(&self) -> f64 {
        -(self.c_a.ln()*self.c_a + self.c_b.ln()*self.c_b)
    }
}

impl RandAtom for BinAtoms {
    type C = BinConcentration;

    fn uniform(rng: &mut Pcg64) -> Self {
        if rng.gen() {
            Self::A
        } else {
            Self::B
        }
    }

    fn with_concentration(rng: &mut Pcg64, c: Self::C) -> Self {
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
