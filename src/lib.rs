#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use std::ops::Index;
use std::ops::IndexMut;

use rand::Rng;
use rand::SeedableRng;
use rand_pcg::Pcg64;
use rand_seeder::Seeder;

mod modular_array;

mod binary;
mod ternary;

pub use modular_array::ModularArray;

pub use binary::{BinAtoms, BinConcentration};
pub use ternary::{TerAtoms, TerConcentration};

/// This trait allows for atoms to be chosen uniformly or after a concentration.
pub trait RandAtom {
    type Concentration: Copy;
    fn uniform(rng: &mut Pcg64) -> Self;
    fn with_concentration(rng: &mut Pcg64, c: Self::Concentration) -> Self;
}

#[derive(Clone)]
pub struct ArrayLatice<A, const WIDTH: usize, const HEIGHT: usize>
where
    [(); WIDTH * HEIGHT]:,
{
    energies: fn(A, A) -> f32,
    grid: ModularArray<A, WIDTH, HEIGHT>,
    rng: Pcg64,
    tot_energy: Option<f32>,
}

/// all constructors
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

/// everything energies
impl<A, const WIDTH: usize, const HEIGHT: usize> ArrayLatice<A, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
    A: Copy + Default + RandAtom,
{
    /// This function returns the total energy of the system.
    /// This is fast when the energy is already calculated and recalculates it if it is not.
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

    /// This function returns the local energy around the idx if it was swapped to atom_at_idx
    fn energies_around(&self, idx: (isize, isize), atom_at_idx: A) -> f32 {
        (self.energies)(atom_at_idx, self.grid[(idx.0 + 1, idx.1)])
            + (self.energies)(atom_at_idx, self.grid[(idx.0 - 1, idx.1)])
            + (self.energies)(atom_at_idx, self.grid[(idx.0, idx.1 + 1)])
            + (self.energies)(atom_at_idx, self.grid[(idx.0, idx.1 - 1)])
    }

    /// This function calculates the energy difference when swapping the indexes
    fn calc_delta_e(&mut self, idx_1: (isize, isize), idx_2: (isize, isize)) -> f32 {
        let e_0 = self.energies_around(idx_1, self.grid[idx_1])
            + self.energies_around(idx_2, self.grid[idx_2]);
        let e_1 = self.energies_around(idx_1, self.grid[idx_2])
            + self.energies_around(idx_2, self.grid[idx_1]);
        e_1 - e_0
    }

    /// This function updates the energy if already calculated and recalculates the whole energy
    /// if it is not already calculated.
    fn update_energy(&mut self, delta_e: f32) {
        match self.tot_energy.as_mut() {
            Some(energy) => *energy += delta_e,
            None => {
                self.tot_energy();
            }
        }
    }
}

/// Choosing and swaping indexes
impl<A, const WIDTH: usize, const HEIGHT: usize> ArrayLatice<A, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
    A: Copy + Default + RandAtom,
{
    /// This function chooses two locations in the grid uniformly.
    fn choose_idxs_uniformly(&mut self) -> ((isize, isize), (isize, isize)) {
        let idx_1 = (
            self.rng.gen_range(0..(WIDTH as isize)),
            self.rng.gen_range(0..(HEIGHT as isize)),
        );
        let idx_2 = (
            self.rng.gen_range(0..(WIDTH as isize)),
            self.rng.gen_range(0..(HEIGHT as isize)),
        );
        (idx_1, idx_2)
    }

    /// This function chooses the first index uniformly and then chooses the second index by sampling the distribution twice.
    /// Once for the x offset the second time for the y offset
    fn choose_idxs_with_distribution<T>(&mut self, distr: T) -> ((isize, isize), (isize, isize))
    where
        T: rand::distributions::Distribution<isize> + Copy,
    {
        let idx_1 = (
            self.rng.gen_range(0..(WIDTH as isize)),
            self.rng.gen_range(0..(HEIGHT as isize)),
        );
        let dx = self.rng.sample(distr);
        let dy = self.rng.sample(distr);
        let idx_2 = (idx_1.0 + dx, idx_1.1 + dy);
        (idx_1, idx_2)
    }

    /// this function swaps the indexes and updates the energies with the provided value
    fn swap_idxs(&mut self, delta_e: f32, idx_1: (isize, isize), idx_2: (isize, isize)) {
        self.update_energy(delta_e);
        let temp = *self.grid.index_mut(idx_1);
        *self.grid.index_mut(idx_1) = *self.grid.index(idx_2);
        *self.grid.index_mut(idx_2) = temp;
    }
}

/// all swapping processes
impl<A, const WIDTH: usize, const HEIGHT: usize> ArrayLatice<A, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
    A: Copy + Default + RandAtom,
{
    /// This uniformly chooses two latice point and swaps the elements if the resulting energy is lower.
    /// If there is no swap it repeats this process until it succeds.
    pub fn swap_uniform(&mut self) {
        loop {
            let (idx_1, idx_2) = self.choose_idxs_uniformly();
            let delta_e = self.calc_delta_e(idx_1, idx_2);
            if delta_e <= 0.0 {
                self.swap_idxs(delta_e, idx_1, idx_2);
                return;
            }
        }
    }

    /// This function chooses one latice point randomly and the second by pulling dy and dx from the distribution twice
    /// and swap the atoms if the resulting energy is lower.
    /// If there is no swap it repeats this process until it succeds.
    pub fn swap_distr<T>(&mut self, distr: T)
    where
        T: rand::distributions::Distribution<isize> + Copy,
    {
        loop {
            let (idx_1, idx_2) = self.choose_idxs_with_distribution(distr);
            let delta_e = self.calc_delta_e(idx_1, idx_2);
            if delta_e <= 0.0 {
                self.swap_idxs(delta_e, idx_1, idx_2);
                return;
            }
        }
    }

    /// This function performs a monte carlo swap with the boltzman factor beta = 1/(k_B * T)
    pub fn monte_carlo_swap(&mut self, temperature: f32) {
        loop {
            let (idx_1, idx_2) = self.choose_idxs_uniformly();
            let delta_e = self.calc_delta_e(idx_1, idx_2);
            if delta_e <= 0.0 || (self.rng.gen::<f32>() > (-temperature * delta_e).exp()) {
                self.swap_idxs(delta_e, idx_1, idx_2);
                return;
            }
        }
    }

    /// This function performs a monte carlo swap with the boltzman factor beta = 1/(k_B * T)
    /// using a distribution for the distance between the two latice sites
    pub fn monte_carlo_swap_distr<T>(&mut self, distr: T, beta: f32)
    where
        T: rand::distributions::Distribution<isize> + Copy,
    {
        loop {
            let (idx_1, idx_2) = self.choose_idxs_with_distribution(distr);
            let delta_e = self.calc_delta_e(idx_1, idx_2);
            if delta_e <= 0.0 || (self.rng.gen::<f32>() > (-beta * delta_e).exp()) {
                self.swap_idxs(delta_e, idx_1, idx_2);
                return;
            }
        }
    }
}

impl<A, const WIDTH: usize, const HEIGHT: usize> ArrayLatice<A, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
    A: Copy + Default + RandAtom,
{
    pub fn as_bytes<'a>(&'a self) -> &'a [u8]
    where
        &'a ModularArray<A, WIDTH, HEIGHT>: Into<&'a [u8]>,
    {
        (&self.grid).into()
    }
}
