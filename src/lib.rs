use std::{
    collections::{btree_map::Entry, BTreeMap, HashMap},
    ops::{Index, IndexMut},
    process::Command,
};

use rand_pcg::Pcg64;

mod array_2d;
pub use array_2d::Array2d;
mod array_3d;
pub use array_3d::Array3d;
mod fast_array;
pub use fast_array::FastArray;

mod atoms;
pub use atoms::{BinAtom, BinConcentration, Energies, Mark, RandAtom};

mod system;
pub use system::System;

pub mod anim;
pub mod logs;

type MyRng = Pcg64;

pub trait Lattice: Index<Self::Index, Output = Self::Atom> + IndexMut<Self::Index> {
    type Atom: Copy + RandAtom;
    type Index: Copy;
    type Neighbors: AsRef<[Self::Index]>;

    fn fill_value(val: Self::Atom) -> Self;
    fn fill_with_fn(func: &mut impl FnMut(Self::Index) -> Self::Atom) -> Self;

    fn all_neighbors(&self) -> HashMap<(Self::Atom, Self::Atom), u32>;
    fn all_neighbors_to(&self, idx: Self::Index) -> Self::Neighbors;

    fn all_idxs(&self) -> Vec<Self::Index>;
    fn tot_sites(&self) -> usize;

    fn as_flat_slice(&self) -> &[Self::Atom];
    fn as_flat_slice_mut(&mut self) -> &mut [Self::Atom];

    fn random_idx(&self, rng: &mut MyRng) -> Self::Index;
    fn choose_idxs_uniformly(&self, rng: &mut MyRng) -> (Self::Index, Self::Index) {
        (self.random_idx(rng), self.random_idx(rng))
    }
    // fn choose_idxs_with_distribution(
    //     &self,
    //     rng: &mut MyRng,
    //     distr: impl Distribution<Self::Index>,
    // ) -> (Self::Index, Self::Index);
    fn reduce_index(&self, idx: Self::Index) -> Self::Index;
    fn swap_vals(&mut self, idx_1: Self::Index, idx_2: Self::Index) {
        let temp = self[idx_1];
        self[idx_1] = self[idx_2];
        self[idx_2] = temp;
    }
}

pub trait GifFrame: Lattice {
    fn get_frame(&self) -> gif::Frame<'_>;
}

// Todo why BTreeMap?
pub struct ClusterDistribution(BTreeMap<u32, u32>);

impl ClusterDistribution {
    pub fn new() -> Self {
        Self(BTreeMap::new())
    }

    pub fn from_map(map: BTreeMap<u32, u32>) -> Self {
        Self(map)
    }

    pub fn ref_map(&self) -> &BTreeMap<u32, u32> {
        &self.0
    }

    pub fn as_map(self) -> BTreeMap<u32, u32> {
        self.0
    }

    pub fn combine(&mut self, other: &Self) {
        for (key, val) in other.0.iter() {
            let original_value = match self.0.entry(*key) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => entry.insert(0),
            };
            *original_value += *val
        }
    }
}

pub trait ClusterCounter: Lattice
where
    <Self as Lattice>::Atom: Mark,
{
    fn count_clusters(&mut self, atom: Self::Atom) -> ClusterDistribution;

    /// # Safety
    /// this function uses mark() of the underlying atom type
    /// which has to be undone using .unmark()
    unsafe fn mark_region_and_get_size(&mut self, idx: Self::Index) -> u32;
}

impl<T> ClusterCounter for T
where
    T: Lattice,
    <T as Lattice>::Atom: Mark,
{
    fn count_clusters(&mut self, atom: Self::Atom) -> ClusterDistribution {
        let mut map = BTreeMap::new();
        for idx in self.all_idxs() {
            if !self[idx].is_marked() && self[idx] == atom {
                // SAFETY:
                // this is safe because we unmark all item in the for_each afterwards
                *map.entry(unsafe { self.mark_region_and_get_size(idx) })
                    .or_insert(0) += 1;
            }
        }
        self.as_flat_slice_mut()
            .iter_mut()
            .for_each(|atom| atom.unmark());
        ClusterDistribution(map)
    }

    /// # Safety
    /// this function uses mark() of the underlying atom type
    /// which has to be undone using .unmark()
    unsafe fn mark_region_and_get_size(&mut self, idx: Self::Index) -> u32 {
        let atom_type = *self[idx];
        let mut stack = vec![idx];
        let mut count = 0;
        while let Some(idx) = stack.pop() {
            if *self[idx] == atom_type {
                // Safety: the unsafe is passed to the caller
                unsafe { self[idx].mark() };
                count += 1;
                stack.extend_from_slice(self.all_neighbors_to(idx).as_ref())
            }
        }
        count
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ClusterStats {
    min: u32,
    quart_1: u32,
    median: u32,
    quart_3: u32,
    max: u32,
}

impl ClusterStats {
    pub fn from_map_block(map: &ClusterDistribution) -> Self {
        let tot_blocks = map
            .0
            .iter()
            .fold(0, |acc_count, (_, block_count)| acc_count + block_count);
        let mut vec: Vec<(u32, u32)> = map.0.iter().map(|(a, b)| (*a, *b)).collect();
        vec.sort_by_key(|(size, _count)| *size);
        let mut count_i = 0;
        let mut quart_1 = 0;
        let mut median = 0;
        let mut quart_3 = 0;
        for (size, count) in &vec {
            count_i += count;
            if count_i <= tot_blocks / 4 {
                quart_1 = *size
            }
            if count_i <= tot_blocks / 2 {
                median = *size
            }
            if count_i <= 3 * tot_blocks / 4 {
                quart_3 = *size
            }
        }
        Self {
            min: vec[0].0,
            quart_1,
            median,
            quart_3,
            max: vec.pop().expect("should not be empty").0,
        }
    }

    pub fn from_map_atom(map: &ClusterDistribution) -> Self {
        let tot_atoms = map
            .0
            .iter()
            .fold(0, |acc_count, (block_size, block_count)| {
                acc_count + block_count * block_size
            });
        let mut vec: Vec<(u32, u32)> = map.0.iter().map(|(a, b)| (*a, *b)).collect();
        vec.sort_by_key(|(size, _count)| *size);
        let mut count_i = 0;
        let mut quart_1 = vec[0].0;
        let mut median = vec[0].0;
        let mut quart_3 = vec[0].0;
        for (size, count) in &vec {
            count_i += count * size;
            if count_i <= tot_atoms / 4 {
                quart_1 = *size
            }
            if count_i <= tot_atoms / 2 {
                median = *size
            }
            if count_i <= 3 * tot_atoms / 4 {
                quart_3 = *size
            }
            if count_i > 3 * tot_atoms / 4 {
                break;
            }
        }
        Self {
            min: vec[0].0,
            quart_1,
            median,
            quart_3,
            max: vec.pop().expect("should not be empty").0,
        }
    }

    pub fn get_categories(prefix: Option<impl ToString>) -> Vec<String> {
        let mut out =
            vec![
                "min".to_owned(),
                "quart_1".to_owned(),
                "median".to_owned(),
                "quart_3".to_owned(),
                "max".to_owned(),
            ];
        if let Some(prefix) = prefix {
            out.iter_mut().for_each(|string| {
                let mut temp = prefix.to_string();
                temp.push_str(string);
                *string = temp;
            })
        }
        out
    }

    pub fn as_vec_f32(&self) -> Vec<f32> {
        vec![
            self.min as f32,
            self.quart_1 as f32,
            self.median as f32,
            self.quart_3 as f32,
            self.max as f32,
        ]
    }

    pub fn valid(&self) -> bool {
        self.min <= self.quart_1
            && self.quart_1 <= self.median
            && self.median <= self.quart_3
            && self.quart_3 <= self.max
    }
}

pub fn disentangle(cluster_stats: Vec<ClusterStats>) -> [Vec<u32>; 5] {
    [
        cluster_stats.iter().map(|stats| stats.min).collect(),
        cluster_stats.iter().map(|stats| stats.quart_1).collect(),
        cluster_stats.iter().map(|stats| stats.median).collect(),
        cluster_stats.iter().map(|stats| stats.quart_3).collect(),
        cluster_stats.iter().map(|stats| stats.max).collect(),
    ]
}

pub fn run_python(script: &str, arg: &str) {
    match Command::new("/opt/homebrew/Caskroom/miniconda/base/envs/datasc/bin/python")
        .arg(script)
        .arg(arg)
        .output()
    {
        Ok(out) => {
            eprintln!("{}", String::from_utf8_lossy(&out.stderr));
            println!("{}", String::from_utf8_lossy(&out.stdout));
            if !out.status.success() {
                eprintln!("failed to run python! {}", out.status);
            }
        }
        Err(err) => eprintln!("failed to run python command {}", err),
    }
}

pub fn flatten<T, const N: usize, const M: usize>(arr: &[[T; N]; M]) -> &[T] {
    // SAFETY: `self.len() * N` cannot overflow because `self` is
    // already in the address space.
    // SAFETY: `[T]` is layout-identical to `[T; N]`
    if std::mem::size_of::<T>() == 0 {
        panic!()
    }
    unsafe { std::slice::from_raw_parts(arr.as_ptr().cast(), N * M) }
}

pub struct StreamingStats {
    count: u32,
    m_k: f64,
    m_k_1: f64,
    v_k: f64,
    v_k_1: f64,
}

impl StreamingStats {
    // variance after https://math.stackexchange.com/questions/20593/calculate-variance-from-a-stream-of-sample-values
    pub fn new() -> Self {
        Self {
            count: 0,
            m_k: 0.0,
            m_k_1: 0.0,
            v_k: 0.0,
            v_k_1: 0.0,
        }
    }
    // _1 is old
    pub fn add_value(&mut self, x_k: f32) {
        let x_k = x_k as f64;
        self.count += 1;
        self.m_k_1 = self.m_k;
        self.v_k_1 = self.v_k;
        self.m_k = self.m_k_1 + (x_k - self.m_k_1) / self.count as f64;
        self.v_k = self.v_k_1 + (x_k - self.m_k_1) * (x_k - self.m_k);
    }

    pub fn avg(&self) -> f32 {
        self.m_k as f32
    }

    pub fn variance(&self) -> f32 {
        (self.v_k / self.count as f64) as f32
    }
}
