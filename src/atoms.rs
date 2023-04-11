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
pub struct Concentration<const N: usize> {
    cs: [f64; N],
}

impl<const N: usize> Concentration<N> {
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

impl<const N: usize> Concentration<N> {
    pub fn uniform() -> Self {
        Self::new([1.0; N])
    }

    pub fn max_entropy(&self) -> f64 {
        -(self.cs.iter().fold(0.0, |acc, x| acc + x.ln() * x))
    }
}

impl<const N: usize> Atom<N> {
    pub fn uniform(rng: &mut rand_pcg::Pcg64) -> Self {
        Self(rng.gen_range(0..(N as u8)))
    }

    pub fn with_concentration(rng: &mut rand_pcg::Pcg64, c: Concentration<N>) -> Self {
        let distr = WeightedIndex::new(&c.cs).unwrap();
        Self(distr.sample(rng) as u8)
    }
}

impl<const W: usize, const H: usize, const N: usize> From<&ModularArray<Atom<N>, W, H>> for &[u8]
where
    [(); W * H]:,
{
    fn from(value: &ModularArray<Atom<N>, W, H>) -> Self {
        unsafe {
            let ptr = value.grid.as_ptr() as *const u8;
            let len = value.grid.len() * std::mem::size_of::<Atom<N>>();
            std::slice::from_raw_parts(ptr, len)
        }
    }
}
