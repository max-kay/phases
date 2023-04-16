#![allow(incomplete_features)]
#![feature(generic_const_exprs, concat_idents)]
use std::ops::Index;
use std::ops::IndexMut;

use rand::distributions::Distribution;
use rand::seq::SliceRandom;
use rand::Rng;
use rand::SeedableRng;
use rand_pcg::Pcg64;
use rand_seeder::Seeder;

mod modular_array;

pub mod anim;
pub mod plots;

mod atoms;
pub use atoms::{AtomNum, ConcentrationNum};

pub use modular_array::ModularArray;

type MyRng = Pcg64;

pub trait Storage: Index<Self::Index, Output = Self::Atom> + IndexMut<Self::Index> {
    type Atom: Copy + RandAtom;
    type Index: Copy;
    fn fill_value(val: Self::Atom) -> Self;
    fn fill_with_fn(func: &mut impl FnMut(Self::Index) -> Self::Atom) -> Self;
    fn all_neighbours(&self) -> Vec<(Self::Atom, Self::Atom)>;
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

pub trait RandAtom: Default + Eq + PartialEq {
    type RandConcentration: Copy + Concentration;
    fn uniform(rng: &mut MyRng) -> Self;
    fn with_concentration(rng: &mut MyRng, cs: Self::RandConcentration) -> Self;
}

pub trait Concentration {
    fn uniform() -> Self;
    fn max_entropy(&self) -> f32;
}

pub struct System<S: Storage> {
    bond_energies: fn(S::Atom, S::Atom) -> f32,
    concentration: <S::Atom as RandAtom>::RandConcentration,
    storage: S,
    rng: MyRng,
    internal_energy: Option<f32>,
    vacancy: Option<S::Index>,
}

/// all constructors
impl<S: Storage> System<S> {
    pub fn new(
        bond_energies: fn(S::Atom, S::Atom) -> f32,
        seed: Option<&str>,
        concentration: Option<<S::Atom as RandAtom>::RandConcentration>,
    ) -> Self {
        let mut rng = match seed {
            Some(seed) => Seeder::from(seed).make_rng(),
            None => MyRng::from_entropy(),
        };

        let (grid, concentration) = match concentration {
            None => (
                S::fill_with_fn(&mut |_| S::Atom::uniform(&mut rng)),
                <S::Atom as RandAtom>::RandConcentration::uniform(),
            ),
            Some(concentration) => (
                S::fill_with_fn(&mut |_| {
                    <S::Atom as RandAtom>::with_concentration(&mut rng, concentration)
                }),
                concentration,
            ),
        };

        let mut obj = Self {
            bond_energies,
            storage: grid,
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
impl<S: Storage> System<S> {
    /// This function returns the total energy of the system.
    /// This is fast when the energy is already calculated and recalculates it if it is not.
    pub fn internal_energy(&mut self) -> f32 {
        if let Some(energy) = self.internal_energy {
            energy
        } else {
            let energy = self
                .storage
                .all_neighbours()
                .iter()
                .fold(0.0, |energy, (a1, a2)| {
                    energy + (self.bond_energies)(*a1, *a2)
                });
            self.internal_energy = Some(energy);
            energy
        }
    }

    /// This function returns the local energy around the idx if it was swapped to atom_at_idx
    fn energies_around(&self, idx: S::Index) -> f32 {
        self.storage
            .all_neighbours_to(idx)
            .iter()
            .fold(0.0, |energy, idx_i| {
                energy + (self.bond_energies)(self.storage[idx], self.storage[*idx_i])
            })
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

impl<S: Storage> System<S> {
    pub fn get_cs(&self) -> <S::Atom as RandAtom>::RandConcentration {
        self.concentration
    }
}

/// all swapping processes
impl<S: Storage> System<S> {
    /// This uniformly chooses two lattice point and swaps the elements if the resulting energy is lower.
    /// If there is no swap it repeats this process until it succeds.
    pub fn swap_uniform(&mut self) -> bool {
        let (idx_1, idx_2) = loop {
            let (idx_1, idx_2) = self.storage.choose_idxs_uniformly(&mut self.rng);
            if self.storage[idx_1] != self.storage[idx_2] {
                break (idx_1, idx_2);
            }
        };
        let e_0 = self.energies_around(idx_1) + self.energies_around(idx_2);
        self.storage.swap_idxs(idx_1, idx_2);
        let e_1 = self.energies_around(idx_1) + self.energies_around(idx_2);
        let delta_e = e_1 - e_0;
        if delta_e <= 0.0 {
            self.update_energy(delta_e);
            true
        } else {
            self.storage.swap_idxs(idx_1, idx_2);
            false
        }
    }

    /// This function chooses one lattice point randomly and the second by pulling dy and dx from the distribution twice
    /// and swap the atoms if the resulting energy is lower.
    /// If there is no swap it repeats this process until it succeds.
    pub fn swap_distr<T>(&mut self, distr: T) -> bool
    where
        T: Distribution<S::Index> + Copy,
    {
        let (idx_1, idx_2) = loop {
            let (idx_1, idx_2) = self
                .storage
                .choose_idxs_with_distribution(&mut self.rng, distr);
            if self.storage[idx_1] != self.storage[idx_2] {
                break (idx_1, idx_2);
            }
        };
        let e_0 = self.energies_around(idx_1) + self.energies_around(idx_2);
        self.storage.swap_idxs(idx_1, idx_2);
        let e_1 = self.energies_around(idx_1) + self.energies_around(idx_2);
        let delta_e = e_1 - e_0;
        if delta_e <= 0.0 {
            self.update_energy(delta_e);
            true
        } else {
            self.storage.swap_idxs(idx_1, idx_2);
            false
        }
    }

    /// This function performs a monte carlo swap with the boltzman factor beta = 1/(k_B * T)
    pub fn monte_carlo_swap(&mut self, beta: f32) -> bool {
        let (idx_1, idx_2) = loop {
            let (idx_1, idx_2) = self.storage.choose_idxs_uniformly(&mut self.rng);
            if self.storage[idx_1] != self.storage[idx_2] {
                break (idx_1, idx_2);
            }
        };
        let e_0 = self.energies_around(idx_1) + self.energies_around(idx_2);
        self.storage.swap_idxs(idx_1, idx_2);
        let e_1 = self.energies_around(idx_1) + self.energies_around(idx_2);
        let delta_e = e_1 - e_0;
        if delta_e <= 0.0 || (self.rng.gen::<f32>() < (-beta * delta_e).exp()) {
            self.update_energy(delta_e);
            true
        } else {
            self.storage.swap_idxs(idx_1, idx_2);
            false
        }
    }

    /// This function performs a monte carlo swap with the boltzman factor beta = 1/(k_B * T)
    /// using a distribution for the distance between the two lattice sites
    pub fn monte_carlo_swap_distr<T>(&mut self, distr: T, beta: f32) -> bool
    where
        T: Distribution<S::Index> + Copy,
    {
        let (idx_1, idx_2) = loop {
            let (idx_1, idx_2) = self
                .storage
                .choose_idxs_with_distribution(&mut self.rng, distr);
            if self.storage[idx_1] != self.storage[idx_2] {
                break (idx_1, idx_2);
            }
        };
        let e_0 = self.energies_around(idx_1) + self.energies_around(idx_2);
        self.storage.swap_idxs(idx_1, idx_2);
        let e_1 = self.energies_around(idx_1) + self.energies_around(idx_2);
        let delta_e = e_1 - e_0;
        if delta_e <= 0.0 || (self.rng.gen::<f32>() < (-beta * delta_e).exp()) {
            self.update_energy(delta_e);
            true
        } else {
            self.storage.swap_idxs(idx_1, idx_2);
            false
        }
    }

    pub fn move_vacancy(&mut self, beta: f32) -> bool {
        // this function does not put any other values in to the vacancy spot.
        // it just reinterprets this spot as being empty thus having bond energies = 0

        if let Some(idx) = self.vacancy {
            let all_neighbours_to = self.storage.all_neighbours_to(idx);
            let other_idx = all_neighbours_to
                .choose(&mut self.rng)
                .expect("all atoms have neighbours");

            // TODO argue why this is correct
            let e_0 = self.energies_around(*other_idx);
            self.storage.swap_idxs(idx, *other_idx);
            let e_1 = self.energies_around(idx);

            let delta_e = e_1 - e_0;
            if delta_e <= 0.0 || (self.rng.gen::<f32>() < (-beta * delta_e).exp()) {
                self.update_energy(delta_e);
                self.vacancy = Some(*other_idx);
                true
            } else {
                self.storage.swap_idxs(idx, *other_idx);
                false
            }
        } else {
            self.vacancy = Some(self.storage.random_idx(&mut self.rng));
            self.move_vacancy(beta)
        }
    }
}

impl<const N: usize, const W: usize, const H: usize> System<ModularArray<AtomNum<N>, W, H>>
where
    [(); W * H]:,
{
    pub fn as_bytes<'a>(&'a self) -> &'a [u8]
    where
        &'a ModularArray<AtomNum<N>, W, H>: Into<&'a [u8]>,
    {
        // how should I deal with vacancies TODO
        (&self.storage).into()
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
