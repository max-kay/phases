use std::{
    collections::HashMap,
    ops::{Index, IndexMut},
};

use rand::Rng;

use crate::{ATrait, Lattice, MyRng, NumAtom, RegionCounter};

/// A 2D grid type that is Copy and allows indexes to "wrap around"
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

impl<T: Copy + ATrait, const W: usize, const H: usize> Lattice for Array2d<T, W, H> {
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

    fn random_idx(&self, rng: &mut MyRng) -> Self::Index {
        (rng.gen_range(0..W as isize), rng.gen_range(0..H as isize))
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

impl<const W: usize, const H: usize, const N: usize> RegionCounter for Array2d<NumAtom<N>, W, H> {
    fn count_regions(&self) -> HashMap<Self::Atom, HashMap<Self::Atom, u32>> {
        todo!()
    }
}

// Box<[[T; W]; H]>
impl<const W: usize, const H: usize, const N: usize> Array2d<NumAtom<N>, W, H> {
    pub fn get_slice(&self) -> &[u8] {
        let ptr: *const NumAtom<N> = &self.grid[0][0];
        let len = H * W * std::mem::size_of::<NumAtom<N>>();
        // SAFETY: This is safe because we know the lenght of the array and NumAtom<N>
        // always just contains an u8 and is repr(transparent)
        // and since we have &self there is nothing mutating the array
        unsafe { std::slice::from_raw_parts(ptr as *const u8, len) }
    }
}
