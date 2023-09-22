use std::{
    collections::HashMap,
    ops::{Index, IndexMut},
};

use rand::Rng;

use crate::{BinAtom, GifFrame, Lattice, RandAtom};

/// A modular 2d Array this implementation uses bit manipulation to implement the modularity of the grid.
/// Because of this SIDE needs to be a power of two with POW beeing this power.
/// (SIDE == 2^POW)
pub struct FastArray<T: Copy, const SIDE: usize, const POW: usize>(Box<[[T; SIDE]; SIDE]>);

impl<T: Copy, const SIDE: usize, const POW: usize> FastArray<T, SIDE, POW> {
    const MASK: usize = !(usize::MAX << POW);

    unsafe fn get_unchecked(&self, idx: (usize, usize)) -> &T {
        self.0.get_unchecked(idx.1).get_unchecked(idx.0)
    }

    unsafe fn get_unchecked_mut(&mut self, idx: (usize, usize)) -> &mut T {
        self.0.get_unchecked_mut(idx.1).get_unchecked_mut(idx.0)
    }
}

impl<T: Copy + Default, const SIDE: usize, const POW: usize> FastArray<T, SIDE, POW> {
    pub fn new() -> Self {
        debug_assert!(SIDE == 2_usize.pow(POW as u32));
        Self(Box::new([[Default::default(); SIDE]; SIDE]))
    }
}

impl<T: Copy, const SIDE: usize, const POW: usize> Index<(usize, usize)>
    for FastArray<T, SIDE, POW>
{
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let x = index.0 & Self::MASK;
        let y = index.1 & Self::MASK;
        debug_assert!(x < SIDE && y < SIDE);
        // Safety: the mask masks all bits that should be zero
        unsafe { self.get_unchecked((x, y)) }
    }
}

impl<T: Copy, const SIDE: usize, const POW: usize> IndexMut<(usize, usize)>
    for FastArray<T, SIDE, POW>
{
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let x = index.0 & Self::MASK;
        let y = index.1 & Self::MASK;
        debug_assert!(x < SIDE && y < SIDE);
        // Safety: the mask masks all bits that should be zero
        unsafe { self.get_unchecked_mut((x, y)) }
    }
}

impl<T: Copy + Default, const SIDE: usize, const POW: usize> Default for FastArray<T, SIDE, POW> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Copy + RandAtom, const SIDE: usize, const POW: usize> Lattice for FastArray<T, SIDE, POW> {
    type Atom = T;

    type Index = (usize, usize);

    type Neighbors = [Self::Index; 4];

    fn fill_value(val: Self::Atom) -> Self {
        Self(Box::new([[val; SIDE]; SIDE]))
    }

    fn fill_with_fn(func: &mut impl FnMut(Self::Index) -> Self::Atom) -> Self {
        let mut new = Self::new();
        for y in 0..SIDE {
            for x in 0..SIDE {
                // Safety: by the above loop x and y are < SIDE
                unsafe { *new.get_unchecked_mut((x, y)) = func((x, y)) }
            }
        }
        new
    }

    fn all_neighbors(&self) -> HashMap<(Self::Atom, Self::Atom), u32> {
        let mut out = HashMap::new();
        for y in 0..SIDE {
            for x in 0..SIDE {
                // Safety: by the above loop x and y are < SIDE
                unsafe {
                    out.entry((*self.get_unchecked((x, y)), self[(x + 1, y)]))
                        .and_modify(|e| *e += 1)
                        .or_insert(1);

                    out.entry((*self.get_unchecked((x, y)), self[(x, y + 1)]))
                        .and_modify(|e| *e += 1)
                        .or_insert(1);
                }
            }
        }
        out
    }

    fn all_neighbors_to(&self, idx: Self::Index) -> Self::Neighbors {
        let x = idx.0;
        let y = idx.1;
        [
            (x.wrapping_add(1), y),
            (x, y.wrapping_add(1)),
            (x.wrapping_add_signed(-1), y),
            (x, y.wrapping_add_signed(-1)),
        ]
    }

    fn all_idxs(&self) -> Vec<Self::Index> {
        let mut out = Vec::with_capacity(SIDE * SIDE);
        for y in 0..SIDE {
            for x in 0..SIDE {
                out.push((x, y));
            }
        }
        out
    }

    fn tot_sites(&self) -> usize {
        SIDE * SIDE
    }

    fn as_flat_slice(&self) -> &[Self::Atom] {
        // Safety: the memory layout of [[T; N]; N] is the same as [T]
        // TODO Zero sized types
        unsafe { std::slice::from_raw_parts(self.0.as_ptr().cast(), SIDE * SIDE) }
    }

    fn as_flat_slice_mut(&mut self) -> &mut [Self::Atom] {
        // Safety: the memory layout of [[T; N]; N] is the same as [T]
        // TODO Zero sized types
        unsafe { std::slice::from_raw_parts_mut(self.0.as_mut_ptr().cast(), SIDE * SIDE) }
    }

    fn random_idx(&self, rng: &mut crate::MyRng) -> Self::Index {
        (rng.gen_range(0..SIDE), rng.gen_range(0..SIDE))
    }

    fn reduce_index(&self, idx: Self::Index) -> Self::Index {
        (idx.0 & Self::MASK, idx.0 & Self::MASK)
    }
}

impl<const SIDE: usize, const POW: usize> GifFrame for FastArray<BinAtom, SIDE, POW> {
    fn get_frame(&self) -> gif::Frame<'_> {
        gif::Frame::from_indexed_pixels(
            SIDE as u16,
            SIDE as u16,
            // Safety: the memory layout of [[T; N]; N] is the same as [T]
            unsafe { std::slice::from_raw_parts(self.0.as_ptr().cast(), SIDE * SIDE) },
            None,
        )
    }
}
