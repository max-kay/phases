#![feature(concat_idents)]
use std::{
    collections::HashMap,
    hash::Hash,
    ops::{Index, IndexMut},
    process::Command,
};

use rand::distributions::Distribution;
use rand_pcg::Pcg64;

mod array_2d;
pub use array_2d::Array2d;
mod array_3d;
pub use array_3d::Array3d;

mod atoms;
pub use atoms::{get_energies_dict, NumAtom, NumC};

mod system;
pub use system::System;

pub mod anim;
pub mod logs;

type MyRng = Pcg64;

pub trait Lattice: Index<Self::Index, Output = Self::Atom> + IndexMut<Self::Index> {
    type Atom: Copy + ATrait;
    type Index: Copy;
    fn fill_value(val: Self::Atom) -> Self;
    fn fill_with_fn(func: &mut impl FnMut(Self::Index) -> Self::Atom) -> Self;

    fn all_neighbours(&self) -> HashMap<(Self::Atom, Self::Atom), u32>;
    fn all_neighbours_to(&self, idx: Self::Index) -> Vec<Self::Index>;

    fn random_idx(&self, rng: &mut MyRng) -> Self::Index;
    fn choose_idxs_uniformly(&self, rng: &mut MyRng) -> (Self::Index, Self::Index) {
        (self.random_idx(rng), self.random_idx(rng))
    }
    fn choose_idxs_with_distribution(
        &self,
        rng: &mut MyRng,
        distr: impl Distribution<Self::Index>,
    ) -> (Self::Index, Self::Index);
    fn reduce_index(&self, idx: Self::Index) -> Self::Index;
    fn swap_idxs(&mut self, idx_1: Self::Index, idx_2: Self::Index);
}

pub trait ATrait: Default + Eq + PartialEq + Hash {
    type Concentration: Copy + CTrait;
    fn uniform(rng: &mut MyRng) -> Self;
    fn with_concentration(rng: &mut MyRng, cs: Self::Concentration) -> Self;
}

pub trait CTrait {
    fn uniform() -> Self;
    fn max_entropy(&self) -> f32;
}

pub fn run_python(script: &str, arg: &str) {
    match Command::new("/opt/homebrew/Caskroom/miniconda/base/envs/datasc/bin/python")
        .arg(script)
        .arg(arg)
        .output()
    {
        Ok(out) => {
            eprintln!("{}", String::from_utf8_lossy(&out.stderr));
            println!("{}", String::from_utf8_lossy(&out.stdout));
            if !out.status.success() {
                eprintln!("failed to run python! {}", out.status);
            }
        }
        Err(err) => eprintln!("failed to run python command {}", err),
    }
}
