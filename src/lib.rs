#![allow(incomplete_features)]
#![feature(generic_const_exprs, concat_idents)]
use std::ops::Index;
use std::ops::IndexMut;

use rand::distributions::Distribution;
use rand::Rng;
use rand::SeedableRng;
use rand_pcg::Pcg64;
use rand_seeder::Seeder;

mod modular_array;

pub mod anim;
pub mod plots;

mod atoms;
pub use atoms::{Atom, Concentration};

pub use modular_array::ModularArray;

pub struct Lattice<const N: usize, const W: usize, const H: usize>
where
    [(); W * H]:,
{
    bond_energies: fn(Atom<N>, Atom<N>) -> f32,
    concentration: Concentration<N>,
    grid: ModularArray<Atom<N>, W, H>,
    rng: Pcg64,
    internal_energy: Option<f32>,
    vacancy: Option<(isize, isize)>,
}

/// all constructors
impl<const N: usize, const W: usize, const H: usize> Lattice<N, W, H>
where
    [(); W * H]:,
{
    pub fn new(
        bond_energies: fn(Atom<N>, Atom<N>) -> f32,
        seed: Option<&str>,
        concentration: Option<Concentration<N>>,
    ) -> Self {
        let mut rng = match seed {
            Some(seed) => Seeder::from(seed).make_rng(),
            None => Pcg64::from_entropy(),
        };

        let mut grid = ModularArray::new();
        let concentration = match concentration {
            Some(concentration) => {
                for x in 0..W {
                    for y in 0..H {
                        grid[(x, y)] = Atom::with_concentration(&mut rng, concentration);
                    }
                }
                concentration
            }
            None => {
                for x in 0..W {
                    for y in 0..H {
                        grid[(x, y)] = Atom::uniform(&mut rng);
                    }
                }
                Concentration::uniform()
            }
        };

        let mut obj = Self {
            bond_energies,
            grid,
            rng,
            concentration,
            internal_energy: None,
            vacancy: None,
        };
        obj.internal_energy();
        obj
    }
}

/// everything energies
impl<const N: usize, const W: usize, const H: usize> Lattice<N, W, H>
where
    [(); W * H]:,
{
    /// This function returns the total energy of the system.
    /// This is fast when the energy is already calculated and recalculates it if it is not.
    pub fn internal_energy(&mut self) -> f32 {
        if let Some(energy) = self.internal_energy {
            energy
        } else {
            let mut energy = 0.0;
            for x in 0..(W as isize) {
                for y in 0..(H as isize) {
                    energy += (self.bond_energies)(self.grid[(x, y)], self.grid[(x - 1, y)]);
                    energy += (self.bond_energies)(self.grid[(x, y)], self.grid[(x, y - 1)]);
                }
            }
            self.internal_energy = Some(energy);
            energy
        }
    }

    /// This function returns the local energy around the idx if it was swapped to atom_at_idx
    fn energies_around(&self, idx: (isize, isize)) -> f32 {
        (self.bond_energies)(self.grid[idx], self.grid[(idx.0 + 1, idx.1)])
            + (self.bond_energies)(self.grid[idx], self.grid[(idx.0 - 1, idx.1)])
            + (self.bond_energies)(self.grid[idx], self.grid[(idx.0, idx.1 + 1)])
            + (self.bond_energies)(self.grid[idx], self.grid[(idx.0, idx.1 - 1)])
    }

    /// This function updates the energy if already calculated and recalculates the whole energy
    /// if it is not already calculated.
    fn update_energy(&mut self, delta_e: f32) {
        match self.internal_energy.as_mut() {
            Some(energy) => *energy += delta_e,
            None => {
                panic!();
            }
        }
    }
}

/// Choosing and swaping indexes
impl<const N: usize, const W: usize, const H: usize> Lattice<N, W, H>
where
    [(); W * H]:,
{
    /// This function chooses two locations in the grid uniformly.
    fn choose_idxs_uniformly(&mut self) -> ((isize, isize), (isize, isize)) {
        let idx_1 = (
            self.rng.gen_range(0..(W as isize)),
            self.rng.gen_range(0..(H as isize)),
        );
        let idx_2 = (
            self.rng.gen_range(0..(W as isize)),
            self.rng.gen_range(0..(H as isize)),
        );
        (idx_1, idx_2)
    }

    /// This function chooses the first index uniformly and then chooses the second index by sampling the distribution twice.
    /// Once for the x offset the second time for the y offset
    fn choose_idxs_with_distribution<T>(&mut self, distr: T) -> ((isize, isize), (isize, isize))
    where
        T: Distribution<isize> + Copy,
    {
        let idx_1 = (
            self.rng.gen_range(0..(W as isize)),
            self.rng.gen_range(0..(H as isize)),
        );
        let dx = self.rng.sample(distr);
        let dy = self.rng.sample(distr);
        let idx_2 = (idx_1.0 + dx, idx_1.1 + dy);
        (idx_1, idx_2)
    }

    /// this function swaps the indexes
    fn swap_idxs(&mut self, idx_1: (isize, isize), idx_2: (isize, isize)) {
        let temp = *self.grid.index_mut(idx_1);
        *self.grid.index_mut(idx_1) = *self.grid.index(idx_2);
        *self.grid.index_mut(idx_2) = temp;
    }
}

impl<const N: usize, const W: usize, const H: usize> Lattice<N, W, H>
where
    [(); W * H]:,
{
    pub fn get_cs(&self) -> Concentration<N> {
        self.concentration
    }
}

/// all swapping processes
impl<const N: usize, const W: usize, const H: usize> Lattice<N, W, H>
where
    [(); W * H]:,
{
    // should I avoid swapping two atoms of the same type? TODO

    /// This uniformly chooses two lattice point and swaps the elements if the resulting energy is lower.
    /// If there is no swap it repeats this process until it succeds.
    pub fn swap_uniform(&mut self) -> bool {
        let (idx_1, idx_2) = loop {
            let (idx_1, idx_2) = self.choose_idxs_uniformly();
            if self.grid[idx_1] != self.grid[idx_2] {
                break (idx_1, idx_2);
            }
        };
        let e_0 = self.energies_around(idx_1) + self.energies_around(idx_2);
        self.swap_idxs(idx_1, idx_2);
        let e_1 = self.energies_around(idx_1) + self.energies_around(idx_2);
        let delta_e = e_1 - e_0;
        if delta_e <= 0.0 {
            self.update_energy(delta_e);
            true
        } else {
            self.swap_idxs(idx_1, idx_2);
            false
        }
    }

    /// This function chooses one lattice point randomly and the second by pulling dy and dx from the distribution twice
    /// and swap the atoms if the resulting energy is lower.
    /// If there is no swap it repeats this process until it succeds.
    pub fn swap_distr<T>(&mut self, distr: T) -> bool
    where
        T: Distribution<isize> + Copy,
    {
        let (idx_1, idx_2) = loop {
            let (idx_1, idx_2) = self.choose_idxs_with_distribution(distr);
            if self.grid[idx_1] != self.grid[idx_2] {
                break (idx_1, idx_2);
            }
        };
        let e_0 = self.energies_around(idx_1) + self.energies_around(idx_2);
        self.swap_idxs(idx_1, idx_2);
        let e_1 = self.energies_around(idx_1) + self.energies_around(idx_2);
        let delta_e = e_1 - e_0;
        if delta_e <= 0.0 {
            self.update_energy(delta_e);
            true
        } else {
            self.swap_idxs(idx_1, idx_2);
            false
        }
    }

    /// This function performs a monte carlo swap with the boltzman factor beta = 1/(k_B * T)
    pub fn monte_carlo_swap(&mut self, beta: f32) -> bool {
        let (idx_1, idx_2) = loop {
            let (idx_1, idx_2) = self.choose_idxs_uniformly();
            if self.grid[idx_1] != self.grid[idx_2] {
                break (idx_1, idx_2);
            }
        };
        let e_0 = self.energies_around(idx_1) + self.energies_around(idx_2);
        self.swap_idxs(idx_1, idx_2);
        let e_1 = self.energies_around(idx_1) + self.energies_around(idx_2);
        let delta_e = e_1 - e_0;
        if delta_e <= 0.0 || (self.rng.gen::<f32>() < (-beta * delta_e).exp()) {
            self.update_energy(delta_e);
            true
        } else {
            self.swap_idxs(idx_1, idx_2);
            false
        }
    }

    /// This function performs a monte carlo swap with the boltzman factor beta = 1/(k_B * T)
    /// using a distribution for the distance between the two lattice sites
    pub fn monte_carlo_swap_distr<T>(&mut self, distr: T, beta: f32) -> bool
    where
        T: Distribution<isize> + Copy,
    {
        let (idx_1, idx_2) = loop {
            let (idx_1, idx_2) = self.choose_idxs_with_distribution(distr);
            if self.grid[idx_1] != self.grid[idx_2] {
                break (idx_1, idx_2);
            }
        };
        let e_0 = self.energies_around(idx_1) + self.energies_around(idx_2);
        self.swap_idxs(idx_1, idx_2);
        let e_1 = self.energies_around(idx_1) + self.energies_around(idx_2);
        let delta_e = e_1 - e_0;
        if delta_e <= 0.0 || (self.rng.gen::<f32>() < (-beta * delta_e).exp()) {
            self.update_energy(delta_e);
            true
        } else {
            self.swap_idxs(idx_1, idx_2);
            false
        }
    }

    pub fn move_vacancy(&mut self, beta: f32) -> bool {
        // this function does not put any other values in to the vacancy spot.
        // it just reinterprets this spot as being empty thus having bond energies = 0

        if let Some(idx) = self.vacancy {
            let other_idx = match self.rng.gen_range(0..4) {
                0 => (idx.0 + 1, idx.1),
                1 => (idx.0 - 1, idx.1),
                2 => (idx.0, idx.1 + 1),
                3 => (idx.0, idx.1 - 1),
                _ => unreachable!(),
            };

            // this part 
            let e_0 = self.energies_around(other_idx);
            self.swap_idxs(idx, other_idx);
            let e_1 = self.energies_around(idx);


            let delta_e = e_1 - e_0;
            if delta_e <= 0.0 || (self.rng.gen::<f32>() < (-beta * delta_e).exp()) {
                self.update_energy(delta_e);
                self.vacancy = Some(other_idx);
                true
            } else {
                self.swap_idxs(idx, other_idx);
                false
            }
        } else {
            self.vacancy = Some((
                self.rng.gen_range(0..W as isize),
                self.rng.gen_range(0..H as isize),
            ));
            self.move_vacancy(beta)
        }
    }
}

impl<const N: usize, const W: usize, const H: usize> Lattice<N, W, H>
where
    [(); W * H]:,
{
    pub fn as_bytes<'a>(&'a self) -> &'a [u8]
    where
        &'a ModularArray<Atom<N>, W, H>: Into<&'a [u8]>,
    {
        // how should I deal with vacancies
        (&self.grid).into()
    }
}

// #[derive(Default)]
// pub struct Model<const N: usize, const W: usize, const H: usize, D: Distribution<isize>>
// where
//     [(); W * H]:,
// {
//     name: String,
//     concentrations: Vec<Concentration<N>>,
//     energies: Option<fn(Atom<N>, Atom<N>) -> f32>,
//     temps: Vec<f32>,
//     distr: Option<D>,
// }

// impl<const N: usize, const W: usize, const H: usize, D: Distribution<isize>> Model<N, W, H, D>
// where
//     [(); W * H]:,
// {
//     // pub fn new() -> Self {
//     //     // Self::default()
//     // }

//     pub fn run() -> Result<(), Box<dyn Error>> {
//         Ok(())
//     }
// }
