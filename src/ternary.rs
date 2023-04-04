use rand::Rng;
use rand_pcg::Pcg64;

use crate::{Concentration, ModularArray, RandAtom};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
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
    c_b: f64,
    c_c: f64,
    c_b_without_a: f64,
}

impl TerConcentration {
    pub fn new(n_a: f64, n_b: f64, n_c: f64) -> Self {
        let n_tot = n_a + n_b + n_c;
        Self {
            c_a: n_a / n_tot,
            c_b: n_b / n_tot,
            c_c: n_c / n_tot,
            c_b_without_a: n_b / (n_b + n_c),
        }
    }
}

impl Concentration for TerConcentration {
    fn uniform() -> Self {
        Self::new(1.0, 1.0, 1.0)
    }

    fn max_entropy(&self) -> f64 {
        -(self.c_a.ln() * self.c_a + self.c_b.ln() * self.c_b + self.c_c.ln() * self.c_c)
    }
}

impl RandAtom for TerAtoms {
    type C = TerConcentration;

    fn uniform(rng: &mut Pcg64) -> Self {
        match rng.gen_range(0..3) {
            0 => Self::A,
            1 => Self::B,
            2 => Self::C,
            _ => unreachable!(),
        }
    }

    fn with_concentration(rng: &mut Pcg64, c: <Self as RandAtom>::C) -> Self {
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
