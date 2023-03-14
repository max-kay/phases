#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use std::ops::Index;
use std::ops::IndexMut;

use rand::Rng;
use rand::SeedableRng;
use rand_pcg::Pcg64;
use rand_seeder::Seeder;

mod binary;

mod modular_array;
mod trinary;

pub use modular_array::ModularArray;

pub use binary::BinAtoms;
pub use trinary::TriAtoms;

pub trait RandAtom {
    type Concentration: Copy;
    fn uniform(rng: &mut Pcg64) -> Self;
    fn with_concentration(rng: &mut Pcg64, concentration: Self::Concentration) -> Self;
}

#[derive(Clone)]
pub struct ArrayLatice<A, const WIDTH: usize, const HEIGHT: usize>
where
    [(); WIDTH * HEIGHT]:,
{
    energies: fn(A, A) -> f32,
    pub grid: ModularArray<A, WIDTH, HEIGHT>,
    rng: Pcg64,
    tot_energy: Option<f32>,
}

impl<A, const WIDTH: usize, const HEIGHT: usize> ArrayLatice<A, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
    A: Copy + Default + RandAtom,
{
    pub fn new(
        energies: fn(A, A) -> f32,
        seed: Option<&str>,
        concentration: Option<<A as RandAtom>::Concentration>,
    ) -> Self {
        let mut rng = match seed {
            Some(seed) => Seeder::from(seed).make_rng(),
            None => Pcg64::from_entropy(),
        };

        let mut grid = ModularArray::new();
        match concentration {
            Some(concentration) => {
                for x in 0..WIDTH {
                    for y in 0..HEIGHT {
                        grid[(x, y)] = <A as RandAtom>::with_concentration(&mut rng, concentration);
                    }
                }
            }
            None => {
                for x in 0..WIDTH {
                    for y in 0..HEIGHT {
                        grid[(x, y)] = <A as RandAtom>::uniform(&mut rng);
                    }
                }
            }
        }

        Self {
            energies,
            grid,
            rng,
            tot_energy: None,
        }
    }

    pub fn new_from_grid(
        grid: ModularArray<A, WIDTH, HEIGHT>,
        energies: fn(A, A) -> f32,
        seed: Option<&str>,
    ) -> Self {
        let rng = if let Some(seed) = seed {
            Seeder::from(seed).make_rng()
        } else {
            Pcg64::from_entropy()
        };
        Self {
            energies,
            grid,
            rng,
            tot_energy: None,
        }
    }
}
impl<A, const WIDTH: usize, const HEIGHT: usize> ArrayLatice<A, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
    A: Copy + Default + RandAtom,
{
    pub fn tot_energy(&mut self) -> f32 {
        if let Some(energy) = self.tot_energy {
            energy
        } else {
            let mut energy = 0.0;
            for x in 0..(WIDTH as isize) {
                for y in 0..(HEIGHT as isize) {
                    energy += (self.energies)(self.grid[(x, y)], self.grid[(x - 1, y)]);
                    energy += (self.energies)(self.grid[(x, y)], self.grid[(x, y - 1)]);
                }
            }
            self.tot_energy = Some(energy);
            energy
        }
    }

    fn energies_around(&self, idx: (isize, isize), atom_at_idx: A) -> f32 {
        (self.energies)(atom_at_idx, self.grid[(idx.0 + 1, idx.1)])
            + (self.energies)(atom_at_idx, self.grid[(idx.0 - 1, idx.1)])
            + (self.energies)(atom_at_idx, self.grid[(idx.0, idx.1 + 1)])
            + (self.energies)(atom_at_idx, self.grid[(idx.0, idx.1 - 1)])
    }

    fn update_energy(&mut self, delta_e: f32) {
        match self.tot_energy.as_mut() {
            Some(energy) => *energy += delta_e,
            None => {
                self.tot_energy();
            }
        }
    }
}

impl<A, const WIDTH: usize, const HEIGHT: usize> ArrayLatice<A, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
    A: Copy + Default + RandAtom,
{
    /// This uniformly chooses two latice point and swaps the elements if the resultant energy is lower.
    /// If there is no swap it repeats this process until it succeds.
    pub fn swap_uniform(&mut self) {
        loop {
            let idx_1 = (
                self.rng.gen_range(0..(WIDTH as isize)),
                self.rng.gen_range(0..(HEIGHT as isize)),
            );
            let idx_2 = (
                self.rng.gen_range(0..(WIDTH as isize)),
                self.rng.gen_range(0..(HEIGHT as isize)),
            );
            let e_0 = self.energies_around(idx_1, self.grid[idx_1])
                + self.energies_around(idx_2, self.grid[idx_2]);
            let e_1 = self.energies_around(idx_1, self.grid[idx_2])
                + self.energies_around(idx_2, self.grid[idx_1]);
            let delta_e = e_1 - e_0;
            if delta_e <= 0.0 {
                // strange behavior when switched to < TODO
                self.update_energy(delta_e);
                let temp = *self.grid.index_mut(idx_1);
                *self.grid.index_mut(idx_1) = *self.grid.index(idx_2);
                *self.grid.index_mut(idx_2) = temp;
                return;
            }
        }
    }

    /// This function chooses one latice point randomly and the second by pulling dy and dx from the distribution twice
    /// and swap the atoms if the resulting energy is lower.
    /// If there is no swap it repeats this process until it succeds.
    pub fn swap_dist_distr<T>(&mut self, distr: T)
    where
        T: rand::distributions::Distribution<isize> + Copy,
    {
        loop {
            let idx_1 = (
                self.rng.gen_range(0..(WIDTH as isize)),
                self.rng.gen_range(0..(HEIGHT as isize)),
            );

            let dx = self.rng.sample(distr);
            let dy = self.rng.sample(distr);

            let idx_2 = (idx_1.0 + dx, idx_1.1 + dy);

            let e_0 = self.energies_around(idx_1, self.grid[idx_1])
                + self.energies_around(idx_2, self.grid[idx_2]);
            let e_1 = self.energies_around(idx_1, self.grid[idx_2])
                + self.energies_around(idx_2, self.grid[idx_1]);
            let delta_e = e_1 - e_0;
            if delta_e <= 0.0 {
                // does this have strange behavior when switched to < too? TODO YES!!!?
                self.update_energy(delta_e);
                let temp = *self.grid.index_mut(idx_1);
                *self.grid.index_mut(idx_1) = *self.grid.index(idx_2);
                *self.grid.index_mut(idx_2) = temp;
                return;
            }
        }
    }
}
