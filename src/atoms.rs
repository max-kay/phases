use std::{ops::Deref, usize};

use crate::{ATrait, Mark};
use rand::{prelude::*, Rng};
use rand_distr::WeightedIndex;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(transparent)]
/// A newtype around u8 that only allows values smaller than N
pub struct NumAtom<const N: usize>(u8);

impl<const N: usize> Mark for NumAtom<N> {
    /// This constant is for the check if the first bit is allways unoccupied
    const OK: usize = 128 - N;
    /// # Safety
    /// the implementation of this function plays with bits and thus does not garantee
    /// that the value is < N
    /// values have to be unmarked after use
    #[allow(clippy::no_effect, path_statements)]
    unsafe fn mark(&mut self) {
        Self::OK;
        self.0 |= 0b1000_0000
    }

    fn unmark(&mut self) {
        self.0 &= 0b0111_1111
    }

    fn is_marked(&self) -> bool {
        (self.0 & 0b1000_0000) == 0b1000_0000
    }
}

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
    /// Constructor for NumAtom<N> failes if val < N
    pub fn new(val: u8) -> Self {
        assert!(val < N as u8);
        Self(val)
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
