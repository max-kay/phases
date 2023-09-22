use std::{
    collections::HashMap,
    ops::{Index, IndexMut},
};

use itertools::Itertools;
use rand::Rng;

use crate::{BinAtom, GifFrame, Lattice, MyRng, RandAtom};

/// A 3D grid type that is Copy and allows indexes to "wrap around"
#[derive(Clone)]
pub struct Array3d<T, const W: usize, const H: usize, const D: usize> {
    pub grid: Box<[[[T; W]; H]; D]>,
}

impl<T, const W: usize, const H: usize, const D: usize> Array3d<T, W, H, D>
where
    T: Copy,
{
    pub fn new() -> Self
    where
        T: Default + Clone,
    {
        Self {
            grid: Box::new([[[T::default(); W]; H]; D]),
        }
    }
}

impl<T, const W: usize, const H: usize, const D: usize> Index<(isize, isize, isize)>
    for Array3d<T, W, H, D>
{
    type Output = T;

    fn index(&self, index: (isize, isize, isize)) -> &Self::Output {
        let (x, y, z) = index;
        let x = x.rem_euclid(W as isize);
        let y = y.rem_euclid(H as isize);
        let z = z.rem_euclid(D as isize);
        // Safety this is safe because of the above rem
        unsafe {
            (*self.grid)
                .get_unchecked(z as usize)
                .get_unchecked(y as usize)
                .get_unchecked(x as usize)
        }
    }
}

impl<T, const W: usize, const H: usize, const D: usize> IndexMut<(isize, isize, isize)>
    for Array3d<T, W, H, D>
{
    fn index_mut(&mut self, index: (isize, isize, isize)) -> &mut Self::Output {
        let (x, y, z) = index;
        let x = x.rem_euclid(W as isize);
        let y = y.rem_euclid(H as isize);
        let z = z.rem_euclid(D as isize);
        // Safety this is safe because of the above rem
        unsafe {
            (*self.grid)
                .get_unchecked_mut(z as usize)
                .get_unchecked_mut(y as usize)
                .get_unchecked_mut(x as usize)
        }
    }
}

impl<T: Default + Copy, const W: usize, const H: usize, const D: usize> Default
    for Array3d<T, W, H, D>
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Copy + RandAtom, const W: usize, const H: usize, const D: usize> Lattice
    for Array3d<T, W, H, D>
{
    type Atom = T;
    type Index = (isize, isize, isize);
    type Neighbors = [Self::Index; 6];

    fn fill_value(val: Self::Atom) -> Self {
        Self {
            grid: Box::new([[[val; W]; H]; D]),
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

    fn all_neighbors(&self) -> HashMap<(Self::Atom, Self::Atom), u32> {
        let mut out = HashMap::new();
        for x in 0..W as isize {
            for y in 0..H as isize {
                for z in 0..D as isize {
                    out.entry((self[(x, y, z)], self[(x - 1, y, z)]))
                        .and_modify(|e| *e += 1)
                        .or_insert(1);
                    out.entry((self[(x, y, z)], self[(x, y - 1, z)]))
                        .and_modify(|e| *e += 1)
                        .or_insert(1);
                    out.entry((self[(x, y, z)], self[(x, y, z - 1)]))
                        .and_modify(|e| *e += 1)
                        .or_insert(1);
                }
            }
        }
        out
    }

    fn all_neighbors_to(&self, idx: Self::Index) -> Self::Neighbors {
        [
            (idx.0 + 1, idx.1, idx.2),
            (idx.0 - 1, idx.1, idx.2),
            (idx.0, idx.1 + 1, idx.2),
            (idx.0, idx.1 - 1, idx.2),
            (idx.0, idx.1, idx.2 + 1),
            (idx.0, idx.1, idx.2 - 1),
        ]
    }

    fn all_idxs(&self) -> Vec<Self::Index> {
        (0..D as isize)
            .cartesian_product(0..H as isize)
            .cartesian_product(0..W as isize)
            .map(|((z, y), x)| (x, y, z))
            .collect_vec()
    }

    fn as_flat_slice(&self) -> &[Self::Atom] {
        // Safety: the memory layout of [[T; N]; N] is the same as [T]
        // TODO zero sized types
        unsafe { std::slice::from_raw_parts(self.grid.as_ptr().cast(), W * H * D) }
    }

    fn as_flat_slice_mut(&mut self) -> &mut [Self::Atom] {
        // Safety: the memory layout of [[T; N]; N] is the same as [T]
        // TODO zero sized types
        unsafe { std::slice::from_raw_parts_mut(self.grid.as_mut_ptr().cast(), W * H * D) }
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

    fn tot_sites(&self) -> usize {
        W * H * D
    }
}

impl<const W: usize, const H: usize, const D: usize> GifFrame for Array3d<BinAtom, W, H, D> {
    fn get_frame(&self) -> gif::Frame<'_> {
        let ptr = self.grid.as_ptr().cast();
        gif::Frame::from_indexed_pixels(
            W as u16,
            H as u16,
            // Safety: the memory layout of [[T; N]; N] is the same as [T]
            unsafe { std::slice::from_raw_parts(ptr, W * H) },
            None,
        )
    }
}
