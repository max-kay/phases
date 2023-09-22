use std::{
    collections::HashMap,
    ops::{Index, IndexMut},
};

use itertools::Itertools;
use rand::Rng;

use crate::{BinAtom, GifFrame, Lattice, MyRng, RandAtom};

/// A 2D grid type that is Copy and allows indexes to "wrap around"
#[derive(Clone)]
pub struct Array2d<T, const W: usize, const H: usize> {
    pub grid: Box<[[T; W]; H]>,
}

impl<T, const W: usize, const H: usize> Array2d<T, W, H>
where
    T: Copy,
{
    pub fn new() -> Self
    where
        T: Default + Clone,
    {
        Self {
            grid: Box::new([[T::default(); W]; H]),
        }
    }
}

impl<T, const W: usize, const H: usize> Index<(isize, isize)> for Array2d<T, W, H> {
    type Output = T;

    fn index(&self, index: (isize, isize)) -> &Self::Output {
        let (x, y) = index;
        let x = x.rem_euclid(W as isize);
        let y = y.rem_euclid(H as isize);
        // Safety this is safe because of the above rem
        unsafe {
            (*self.grid)
                .get_unchecked(y as usize)
                .get_unchecked(x as usize)
        }
    }
}

impl<T, const W: usize, const H: usize> IndexMut<(isize, isize)> for Array2d<T, W, H> {
    fn index_mut(&mut self, index: (isize, isize)) -> &mut Self::Output {
        let (x, y) = index;
        let x = x.rem_euclid(W as isize);
        let y = y.rem_euclid(H as isize);
        // Safety this is safe because of the above rem
        unsafe {
            (*self.grid)
                .get_unchecked_mut(y as usize)
                .get_unchecked_mut(x as usize)
        }
    }
}

impl<T: Default + Copy, const W: usize, const H: usize> Default for Array2d<T, W, H> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Copy + RandAtom, const W: usize, const H: usize> Lattice for Array2d<T, W, H> {
    type Atom = T;
    type Index = (isize, isize);
    type Neighbors = [Self::Index; 4];

    fn fill_value(val: Self::Atom) -> Self {
        Self {
            grid: Box::new([[val; W]; H]),
        }
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

    fn all_neighbors(&self) -> HashMap<(Self::Atom, Self::Atom), u32> {
        let mut out = HashMap::new();
        for x in 0..W as isize {
            for y in 0..H as isize {
                out.entry((self[(x, y)], self[(x - 1, y)]))
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
                out.entry((self[(x, y)], self[(x, y - 1)]))
                    .and_modify(|e| *e += 1)
                    .or_insert(1);
            }
        }
        out
    }

    fn all_neighbors_to(&self, idx: Self::Index) -> Self::Neighbors {
        [
            (idx.0 + 1, idx.1),
            (idx.0 - 1, idx.1),
            (idx.0, idx.1 + 1),
            (idx.0, idx.1 - 1),
        ]
    }

    fn all_idxs(&self) -> Vec<Self::Index> {
        (0..H as isize)
            .cartesian_product(0..W as isize)
            .map(|(y, x)| (x, y))
            .collect_vec()
    }

    fn as_flat_slice(&self) -> &[Self::Atom] {
        // Safety: the memory layout of [[T; N]; N] is the same as [T]
        // TODO Zero sized types
        unsafe { std::slice::from_raw_parts(self.grid.as_ptr().cast(), W * H) }
    }

    fn as_flat_slice_mut(&mut self) -> &mut [Self::Atom] {
        // Safety: the memory layout of [[T; N]; N] is the same as [T]
        // TODO Zero sized types
        unsafe { std::slice::from_raw_parts_mut(self.grid.as_mut_ptr().cast(), W * H) }
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

    fn tot_sites(&self) -> usize {
        W * H
    }
}

impl<const W: usize, const H: usize> GifFrame for Array2d<BinAtom, W, H> {
    fn get_frame(&self) -> gif::Frame<'_> {
        gif::Frame::from_indexed_pixels(
            W as u16,
            H as u16,
            // Safety: the memory layout of [[T; N]; N] is the same as [T]
            unsafe { std::slice::from_raw_parts(self.grid.as_ptr().cast(), W * H) },
            None,
        )
    }
}
