use std::ops::{Index, IndexMut};

use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;
use rand_seeder::Seeder;

use crate::{ModularArray, RandAtom};

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
    pub fn swap_rand(&mut self) {
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
                self.update_energy(delta_e);
                let temp = *self.grid.index_mut(idx_1);
                *self.grid.index_mut(idx_1) = *self.grid.index(idx_2);
                *self.grid.index_mut(idx_2) = temp;
                return;
            }
        }
    }

    fn energies_around(&self, idx: (isize, isize), atom_at_idx: A) -> f32 {
        (self.energies)(atom_at_idx, self.grid[(idx.0 + 1, idx.1)])
            + (self.energies)(atom_at_idx, self.grid[(idx.0 - 1, idx.1)])
            + (self.energies)(atom_at_idx, self.grid[(idx.0, idx.1 + 1)])
            + (self.energies)(atom_at_idx, self.grid[(idx.0, idx.1 - 1)])
    }
}
