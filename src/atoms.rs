use core::panic;
use std::ops::Deref;

use crate::{ATrait, CTrait};
use rand::{prelude::*, Rng};
use rand_distr::WeightedIndex;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct NumAtom<const N: usize>(u8);

impl<const N: usize> Deref for NumAtom<N> {
    type Target = u8;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const N: usize> ATrait for NumAtom<N> {
    type Concentration = NumC<N>;

    fn uniform(rng: &mut crate::MyRng) -> Self {
        Self(rng.gen_range(0..N as u8))
    }

    fn with_concentration(rng: &mut crate::MyRng, cs: Self::Concentration) -> Self {
        let distr = WeightedIndex::new(&cs.cs).unwrap();
        Self(distr.sample(rng) as u8)
    }
}

impl<const N: usize> NumAtom<N> {
    pub fn new(val: u8) -> Self {
        if val < N as u8 {
            Self(val)
        } else {
            panic!()
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct NumC<const N: usize> {
    cs: [f64; N],
}

impl<const N: usize> NumC<N> {
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

impl<const N: usize> CTrait for NumC<N> {
    fn uniform() -> Self {
        Self::new([1.0; N])
    }
}

pub fn get_energies_dict<const N: usize>(e_func: fn(NumAtom<N>, NumAtom<N>) -> f32) -> String {
    let mut string = "{".to_owned();
    for i in 0..N as u8 {
        for j in i..N as u8 {
            string.push_str(&format!(
                "({}, {}): {}, ",
                i,
                j,
                (e_func)(NumAtom::<N>::new(i), NumAtom::<N>::new(j))
            ))
        }
    }
    string.push('}');
    string
}
