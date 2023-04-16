use core::panic;
use std::ops::Deref;

use crate::Concentration;
use crate::ModularArray;
use crate::RandAtom;
use rand::prelude::*;
use rand::Rng;
use rand_distr::WeightedIndex;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[repr(transparent)]
pub struct AtomNum<const N: usize>(u8);

impl<const N: usize> Deref for AtomNum<N> {
    type Target = u8;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const N: usize> RandAtom for AtomNum<N> {
    type RandConcentration = ConcentrationNum<N>;

    fn uniform(rng: &mut crate::MyRng) -> Self {
        Self(rng.gen_range(0..N as u8))
    }

    fn with_concentration(rng: &mut crate::MyRng, cs: Self::RandConcentration) -> Self {
        let distr = WeightedIndex::new(&cs.cs).unwrap();
        Self(distr.sample(rng) as u8)
    }
}

impl<const N: usize> AtomNum<N> {
    pub fn new(val: u8) -> Self {
        if val < N as u8 {
            Self(val)
        } else {
            panic!()
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ConcentrationNum<const N: usize> {
    cs: [f64; N],
}

impl<const N: usize> ConcentrationNum<N> {
    pub fn new(ns: [f64; N]) -> Self {
        let n_tot: f64 = ns.iter().sum();
        let mut cs = ns;
        cs.iter_mut().for_each(|n| *n /= n_tot);
        Self { cs }
    }

    pub fn get_cs(&self) -> &[f64; N] {
        &self.cs
    }
}

impl<const N: usize> Concentration for ConcentrationNum<N> {
    fn uniform() -> Self {
        Self::new([1.0; N])
    }
    fn max_entropy(&self) -> f32 {
        -(self.cs.iter().fold(0.0, |acc, x| acc + x.ln() * x)) as f32
    }
}

impl<const W: usize, const H: usize, const N: usize> From<&ModularArray<AtomNum<N>, W, H>>
    for &[u8]
where
    [(); W * H]:,
{
    fn from(value: &ModularArray<AtomNum<N>, W, H>) -> Self {
        unsafe {
            let ptr = value.grid.as_ptr() as *const u8;
            let len = value.grid.len() * std::mem::size_of::<AtomNum<N>>();
            std::slice::from_raw_parts(ptr, len)
        }
    }
}
