use std::ops::{Index, IndexMut};


/// A 2D grid type that is Copy and allows indexes to "wrap around" if they're isize
/// and directly acceses the underlying array when using usize
#[derive(Clone, Copy)]
pub struct ModularArray<T, const W: usize, const H: usize>
where
    [(); W * H]:,
{
    pub grid: [T; W * H],
}

impl<T, const W: usize, const H: usize> ModularArray<T, W, H>
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

impl<T, const W: usize, const H: usize> Index<(isize, isize)>
    for ModularArray<T, W, H>
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

impl<T, const W: usize, const H: usize> IndexMut<(isize, isize)>
    for ModularArray<T, W, H>
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

impl<T, const W: usize, const H: usize> Index<(usize, usize)>
    for ModularArray<T, W, H>
where
    [(); W * H]:,
{
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (x, y) = index;
        self.grid.index(y * W + x)
    }
}

impl<T, const W: usize, const H: usize> IndexMut<(usize, usize)>
    for ModularArray<T, W, H>
where
    [(); W * H]:,
{
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (x, y) = index;
        self.grid.index_mut(y * W + x)
    }
}

impl<T: Default + Copy, const W: usize, const H: usize> Default
    for ModularArray<T, W, H>
where
    [(); W * H]:,
{
    fn default() -> Self {
        Self::new()
    }
}
