use core::panic;
use std::ops::Deref;

use crate::{ATrait, Array2d, CTrait};
use rand::{prelude::*, Rng};
use rand_distr::WeightedIndex;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(transparent)]
pub struct NumAtom<const N: u8>(u8);

impl<const N: u8> Deref for NumAtom<N> {
    type Target = u8;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const N: u8> ATrait for NumAtom<N>
where
    [(); N as usize]:,
{
    type Concentration = NumC<N>;

    fn uniform(rng: &mut crate::MyRng) -> Self {
        Self(rng.gen_range(0..N))
    }

    fn with_concentration(rng: &mut crate::MyRng, cs: Self::Concentration) -> Self {
        let distr = WeightedIndex::new(&cs.cs).unwrap();
        Self(distr.sample(rng) as u8)
    }
}

impl<const N: u8> NumAtom<N> {
    pub fn new(val: u8) -> Self {
        if val < N {
            Self(val)
        } else {
            panic!()
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct NumC<const N: u8>
where
    [(); N as usize]:,
{
    cs: [f64; N as usize],
}

impl<const N: u8> NumC<N>
where
    [(); N as usize]:,
{
    pub fn new(ns: [f64; N as usize]) -> Self {
        let n_tot: f64 = ns.iter().sum();
        let mut cs = ns;
        cs.iter_mut().for_each(|n| *n /= n_tot);
        Self { cs }
    }

    pub fn get_cs(&self) -> &[f64; N as usize] {
        &self.cs
    }
}

impl<const N: u8> CTrait for NumC<N>
where
    [(); N as usize]:,
{
    fn uniform() -> Self {
        Self::new([1.0; N as usize])
    }
    fn max_entropy(&self) -> f32 {
        -(self.cs.iter().fold(0.0, |acc, x| acc + x.ln() * x)) as f32
    }
}

impl<const W: usize, const H: usize, const N: u8> From<&Array2d<NumAtom<N>, W, H>> for &[u8]
where
    [(); W * H]:,
{
    fn from(value: &Array2d<NumAtom<N>, W, H>) -> Self {
        unsafe {
            let ptr = value.grid.as_ptr() as *const u8;
            let len = value.grid.len() * std::mem::size_of::<NumAtom<N>>();
            std::slice::from_raw_parts(ptr, len)
        }
    }
}

pub fn get_energies_dict<const N: u8>(e_func: fn(NumAtom<N>, NumAtom<N>) -> f32) -> String {
    let mut string = "{".to_owned();
    for i in 0..N {
        for j in i..N {
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
