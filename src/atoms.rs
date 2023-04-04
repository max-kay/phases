// macro_rules! replace_ty {
//     ($_t:ident $sub:ty) => {
//         $sub
//     };
// }

// macro_rules! replace_expr {
//     ($_t:ident $sub:expr) => {
//         $sub
//     };
// }

// macro_rules! count_tts {
//     ($($tts:tt)*) => {0u8 $(+ replace_expr!($tts 1u8))*};
// }

// macro_rules! atom {
//     ($enum_name:ident, $($variant_name:ident),+) => {
//         mod $enum_name {
//             #[derive(Clone, Copy, Default, PartialEq, Eq)]
//             #[repr(u8)]
//             pub enum Atom {
//                 #[default]
//                 $($variant_name),+
//             }
//             pub struct Concentration;
//             impl $crate::RandAtom for Atom {
//                 type Concentration = ($(replace_ty!($variant_name f64)),+);
//                 fn uniform(rng: &mut rand_pcg::Pcg64) -> Self{
//                     use rand::Rng;
//                     unsafe {
//                         *((rng.gen_range(0..count_tts!($($variant_name)+)) as *const u8) as *const Self)
//                     }
//                 }
//                 fn with_concentration(rng: &mut rand_pcg::Pcg64, c: Self::Concentration) -> Self{
//                     todo!()
//                 }

//             }
//         }
//     };
// }

// atom! {bin, ba, ko}

use std::ops::Deref;

use crate::ModularArray;
use rand::prelude::*;
use rand::Rng;
use rand_distr::WeightedIndex;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[repr(transparent)]
pub struct Atom<const N: usize>(u8);

impl<const N: usize> Deref for Atom<N> {
    type Target = u8;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Concrete<const N: usize> {
    cs: [f64; N],
}

impl<const N: usize> Concrete<N> {
    pub fn new(ns: [f64; N]) -> Self {
        let n_tot: f64 = ns.iter().sum();
        let mut cs = ns;
        cs.iter_mut().for_each(|n| *n /= n_tot);
        Self { cs }
    }
}

impl<const N: usize>  Concrete<N> {
    pub fn uniform() -> Self {
        Self::new([1.0; N])
    }

    pub fn max_entropy(&self) -> f64 {
        -(self.cs.iter().fold(0.0, |acc, x| acc + x.ln() * x))
    }
}

impl<const N_ATOMS: usize>  Atom<N_ATOMS> {
    pub fn uniform(rng: &mut rand_pcg::Pcg64) -> Self {
        Self(rng.gen_range(0..(N_ATOMS as u8)))
    }

    pub fn with_concentration(rng: &mut rand_pcg::Pcg64, c: Concrete<N_ATOMS>) -> Self {
        let distr = WeightedIndex::new(&c.cs).unwrap();
        Self(distr.sample(rng) as u8)
    }
}

impl<const WIDTH: usize, const HEIGHT: usize, const N: usize>
    From<&ModularArray<Atom<N>, WIDTH, HEIGHT>> for &[u8]
where
    [(); WIDTH * HEIGHT]:,
{
    fn from(value: &ModularArray<Atom<N>, WIDTH, HEIGHT>) -> Self {
        unsafe {
            let ptr = value.grid.as_ptr() as *const u8;
            let len = value.grid.len() * std::mem::size_of::<Atom<N>>();
            std::slice::from_raw_parts(ptr, len)
        }
    }
}
