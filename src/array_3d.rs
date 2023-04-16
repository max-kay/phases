use std::ops::{Index, IndexMut};

use rand::Rng;

use crate::{Atom, Lattice, MyRng};

/// A 2D grid type that is Copy and allows indexes to "wrap around" if they're isize
/// and directly acceses the underlying array when using usize
#[derive(Clone, Copy)]
pub struct Array3d<T, const W: usize, const H: usize, const D: usize>
where
    [(); W * H * D]:,
{
    pub grid: [T; W * H * D],
}

impl<T, const W: usize, const H: usize, const D: usize> Array3d<T, W, H, D>
where
    [(); W * H * D]:,
    T: Copy,
{
    pub fn new() -> Self
    where
        T: Default + Clone,
    {
        Self {
            grid: [T::default(); W * H * D],
        }
    }
}

impl<T, const W: usize, const H: usize, const D: usize> Index<(isize, isize, isize)>
    for Array3d<T, W, H, D>
where
    [(); W * H * D]:,
{
    type Output = T;

    fn index(&self, index: (isize, isize, isize)) -> &Self::Output {
        let (x, y, z) = index;
        let x = x.rem_euclid(W as isize);
        let y = y.rem_euclid(H as isize);
        let z = z.rem_euclid(D as isize);
        self.grid
            .index(z as usize * D * W + y as usize * W + x as usize)
    }
}

impl<T, const W: usize, const H: usize, const D: usize> IndexMut<(isize, isize, isize)>
    for Array3d<T, W, H, D>
where
    [(); W * H * D]:,
{
    fn index_mut(&mut self, index: (isize, isize, isize)) -> &mut Self::Output {
        let (x, y, z) = index;
        let x = x.rem_euclid(W as isize);
        let y = y.rem_euclid(H as isize);
        let z = z.rem_euclid(D as isize);
        self.grid
            .index_mut(z as usize * D * W + y as usize * W + x as usize)
    }
}

impl<T: Default + Copy, const W: usize, const H: usize, const D: usize> Default
    for Array3d<T, W, H, D>
where
    [(); W * H * D]:,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Copy + Atom, const W: usize, const H: usize, const D: usize> Lattice for Array3d<T, W, H, D>
where
    [(); W * H * D]:,
{
    type Atom = T;

    type Index = (isize, isize, isize);

    fn fill_value(val: Self::Atom) -> Self {
        Self {
            grid: [val; W * H * D],
        }
    }

    fn fill_with_fn(func: &mut impl FnMut(Self::Index) -> Self::Atom) -> Self {
        let mut out = Self::new();
        for x in 0..W as isize {
            for y in 0..H as isize {
                for z in 0..D as isize {
                    out[(x, y, z)] = func((x, y, z));
                }
            }
        }
        out
    }

    fn all_neighbours(&self) -> Vec<(Self::Atom, Self::Atom)> {
        let mut out = Vec::with_capacity(W * H * 2);
        for x in 0..W as isize {
            for y in 0..H as isize {
                for z in 0..D as isize {
                    out.push((self[(x, y, z)], self[(x - 1, y, z)]));
                    out.push((self[(x, y, z)], self[(x, y - 1, z)]));
                    out.push((self[(x, y, z)], self[(x, y, z - 1)]));
                }
            }
        }
        out
    }

    fn all_neighbours_to(&self, idx: Self::Index) -> Vec<Self::Index> {
        vec![
            (idx.0 + 1, idx.1, idx.2),
            (idx.0 - 1, idx.1, idx.2),
            (idx.0, idx.1 + 1, idx.2),
            (idx.0, idx.1 - 1, idx.2),
            (idx.0, idx.1, idx.2 + 1),
            (idx.0, idx.1, idx.2 - 1),
        ]
    }

    fn choose_idxs_with_distribution(
        &self,
        rng: &mut crate::MyRng,
        distr: impl rand_distr::Distribution<Self::Index>,
    ) -> (Self::Index, Self::Index) {
        let idx_1 = self.random_idx(rng);
        let offset = distr.sample(rng);
        (
            idx_1,
            (idx_1.0 + offset.0, idx_1.1 + offset.1, idx_1.2 + offset.2),
        )
    }

    fn random_idx(&self, rng: &mut MyRng) -> Self::Index {
        (
            rng.gen_range(0..W as isize),
            rng.gen_range(0..H as isize),
            rng.gen_range(0..D as isize),
        )
    }

    fn reduce_index(&self, idx: Self::Index) -> Self::Index {
        let (x, y, z) = idx;
        let x = x.rem_euclid(W as isize);
        let y = y.rem_euclid(H as isize);
        let z = z.rem_euclid(D as isize);
        (x, y, z)
    }

    fn swap_idxs(&mut self, idx_1: Self::Index, idx_2: Self::Index) {
        let temp = self[idx_1];
        self[idx_1] = self[idx_2];
        self[idx_2] = temp;
    }
}
