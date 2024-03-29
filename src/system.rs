use rand::{seq::SliceRandom, Rng, SeedableRng};
use rand_seeder::Seeder;

use crate::{
    ClusterCounter, ClusterDistribution, Energies, GifFrame, Lattice, Mark, MyRng, RandAtom,
};

pub struct System<L: Lattice, E: Energies<L::Atom>> {
    bond_energies: E,
    lattice: L,
    rng: MyRng,
    internal_energy: Option<f32>,
    vacancy: Option<L::Index>,
}

/// all constructors
impl<L: Lattice, E: Energies<L::Atom>> System<L, E> {
    pub fn new(
        bond_energies: E,
        seed: Option<&str>,
        concentration: <L::Atom as RandAtom>::Concentration,
    ) -> Self {
        let mut rng =
            match seed {
                Some(seed) => Seeder::from(seed).make_rng(),
                None => MyRng::from_entropy(),
            };

        let grid = L::fill_with_fn(
            &mut |_| <L::Atom as RandAtom>::with_concentration(&mut rng, concentration)
        );

        let mut obj = Self {
            bond_energies,
            lattice: grid,
            rng,
            internal_energy: None,
            vacancy: None,
        };
        obj.internal_energy();
        obj
    }

    pub fn tot_sites(&self) -> usize {
        self.lattice.tot_sites()
    }

    pub fn get_energies_dict(&self) -> String {
        self.bond_energies.as_dict()
    }
}

/// everything energies
impl<L: Lattice, E: Energies<L::Atom>> System<L, E> {
    /// This function returns the total energy of the system.
    /// This is fast when the energy is already calculated and recalculates it if it is not.
    pub fn internal_energy(&mut self) -> f32 {
        if let Some(energy) = self.internal_energy {
            energy
        } else {
            let energy =
                self.lattice
                    .all_neighbors()
                    .iter()
                    .fold(0.0, |energy, ((a1, a2), count)| {
                        energy + self.bond_energies.get_interaction_energy(*a1, *a2) * *count as f32
                    });
            self.internal_energy = Some(energy);
            energy
        }
    }

    /// This function returns the local energy around the idx if it was swapped to atom_at_idx
    fn energies_around(&self, idx: L::Index) -> f32 {
        self.lattice
            .all_neighbors_to(idx)
            .as_ref()
            .iter()
            .fold(0.0, |energy, idx_i| {
                energy
                    + self
                        .bond_energies
                        .get_interaction_energy(self.lattice[idx], self.lattice[*idx_i])
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

/// all swapping processes
impl<L: Lattice, E: Energies<L::Atom>> System<L, E> {
    /// This function performs a monte carlo swap with the boltzman factor beta = 1/(k_B * T)
    pub fn monte_carlo_swap(&mut self, beta: f32) -> bool {
        let (idx_1, idx_2) = loop {
            let (idx_1, idx_2) = self.lattice.choose_idxs_uniformly(&mut self.rng);
            if self.lattice[idx_1] != self.lattice[idx_2] {
                break (idx_1, idx_2);
            }
        };
        let e_0 = self.energies_around(idx_1) + self.energies_around(idx_2);
        self.lattice.swap_vals(idx_1, idx_2);
        let e_1 = self.energies_around(idx_1) + self.energies_around(idx_2);
        let delta_e = e_1 - e_0;
        if delta_e <= 0.0 || (self.rng.gen::<f32>() < (-beta * delta_e).exp()) {
            self.update_energy(delta_e);
            true
        } else {
            self.lattice.swap_vals(idx_1, idx_2);
            false
        }
    }

    /// this function does not put any other values in to the vacancy spot.
    /// it just reinterprets this spot as being empty thus having bond energies = 0
    pub fn move_vacancy(&mut self, beta: f32) -> bool {
        if let Some(idx) = self.vacancy {
            let all_neighbors_to = self.lattice.all_neighbors_to(idx);
            let other_idx = all_neighbors_to
                .as_ref()
                .choose(&mut self.rng)
                .expect("all atoms have neighbors");

            // This is correct because this model assumes there is no interaction between
            // a neighboring atom and a vacancy
            // The value in the grid where the atom sits leads to an over counting
            // of the energy between index 1 and 2
            // but this doesnt matter because it is in e_0 and e_1 and thus subtrackted out
            // !!SEE COMMENT ON Energies IMPLEMENTATION FOR [f32; 4]
            let e_0 = self.energies_around(*other_idx);
            self.lattice.swap_vals(idx, *other_idx);
            let e_1 = self.energies_around(idx);

            let delta_e = e_1 - e_0;
            if delta_e <= 0.0 || (self.rng.gen::<f32>() < (-beta * delta_e).exp()) {
                self.update_energy(delta_e);
                self.vacancy = Some(*other_idx);
                true
            } else {
                self.lattice.swap_vals(idx, *other_idx);
                false
            }
        } else {
            let idx = self.lattice.random_idx(&mut self.rng);
            self.vacancy = Some(idx);
            self.lattice[idx] = L::Atom::vacancy();
            self.move_vacancy(beta)
        }
    }
}

impl<L: GifFrame, E: Energies<L::Atom>> System<L, E> {
    pub fn get_frame(&self) -> gif::Frame<'_> {
        self.lattice.get_frame()
    }
}

impl<L: ClusterCounter, E: Energies<L::Atom>> System<L, E>
where
    <L as Lattice>::Atom: Mark,
{
    pub fn count_all_clusters(&mut self) -> Vec<ClusterDistribution> {
        let mut out = Vec::new();
        for atom in L::Atom::all_atoms() {
            out.push(self.lattice.count_clusters(atom))
        }
        out
    }
    pub fn count_clusters(&mut self, atom: L::Atom) -> ClusterDistribution {
        self.lattice.count_clusters(atom)
    }
}
