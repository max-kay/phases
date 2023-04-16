use std::ops::{Index, IndexMut};

use rand::Rng;

use crate::{Atom, Lattice, MyRng};

/// A 2D grid type that is Copy and allows indexes to "wrap around" if they're isize
/// and directly acceses the underlying array when using usize
#[derive(Clone, Copy)]
pub struct Array2d<T, const W: usize, const H: usize>
where
    [(); W * H]:,
{
    pub grid: [T; W * H],
}

impl<T, const W: usize, const H: usize> Array2d<T, W, H>
where
    [(); W * H]:,
    T: Copy,
{
    pub fn new() -> Self
    where
        T: Default + Clone,
    {
        Self {
            grid: [T::default(); W * H],
        }
    }
}

impl<T, const W: usize, const H: usize> Index<(isize, isize)> for Array2d<T, W, H>
where
    [(); W * H]:,
{
    type Output = T;

    fn index(&self, index: (isize, isize)) -> &Self::Output {
        let (x, y) = index;
        let x = x.rem_euclid(W as isize);
        let y = y.rem_euclid(H as isize);
        self.grid.index(y as usize * W + x as usize)
    }
}

impl<T, const W: usize, const H: usize> IndexMut<(isize, isize)> for Array2d<T, W, H>
where
    [(); W * H]:,
{
    fn index_mut(&mut self, index: (isize, isize)) -> &mut Self::Output {
        let (x, y) = index;
        let x = x.rem_euclid(W as isize);
        let y = y.rem_euclid(H as isize);
        self.grid.index_mut(y as usize * W + x as usize)
    }
}

impl<T: Default + Copy, const W: usize, const H: usize> Default for Array2d<T, W, H>
where
    [(); W * H]:,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Copy + Atom, const W: usize, const H: usize> Lattice for Array2d<T, W, H>
where
    [(); W * H]:,
{
    type Atom = T;

    type Index = (isize, isize);

    fn fill_value(val: Self::Atom) -> Self {
        Self { grid: [val; W * H] }
    }

    fn fill_with_fn(func: &mut impl FnMut(Self::Index) -> Self::Atom) -> Self {
        let mut out = Self::new();
        for x in 0..W as isize {
            for y in 0..H as isize {
                out[(x, y)] = func((x, y));
            }
        }
        out
    }

    fn all_neighbours(&self) -> Vec<(Self::Atom, Self::Atom)> {
        let mut out = Vec::with_capacity(W * H * 2);
        for x in 0..W as isize {
            for y in 0..H as isize {
                out.push((self[(x, y)], self[(x - 1, y)]));
                out.push((self[(x, y)], self[(x, y - 1)]));
            }
        }
        out
    }

    fn all_neighbours_to(&self, idx: Self::Index) -> Vec<Self::Index> {
        vec![
            (idx.0 + 1, idx.1),
            (idx.0 - 1, idx.1),
            (idx.0, idx.1 + 1),
            (idx.0, idx.1 - 1),
        ]
    }

    fn choose_idxs_with_distribution(
        &self,
        rng: &mut crate::MyRng,
        distr: impl rand_distr::Distribution<Self::Index>,
    ) -> (Self::Index, Self::Index) {
        let idx_1 = self.random_idx(rng);
        let offset = distr.sample(rng);
        (idx_1, (idx_1.0 + offset.0, idx_1.1 + offset.1))
    }

    fn random_idx(&self, rng: &mut MyRng) -> Self::Index {
        (rng.gen_range(0..W as isize), rng.gen_range(0..H as isize))
    }

    fn reduce_index(&self, idx: Self::Index) -> Self::Index {
        let (x, y) = idx;
        let x = x.rem_euclid(W as isize);
        let y = y.rem_euclid(H as isize);
        (x, y)
    }

    fn swap_idxs(&mut self, idx_1: Self::Index, idx_2: Self::Index) {
        let temp = self[idx_1];
        self[idx_1] = self[idx_2];
        self[idx_2] = temp;
    }
}
