use std::ops::{Index, IndexMut};

#[derive(Clone, Copy)]
pub struct ModularArray<T, const WIDTH: usize, const HEIGHT: usize>
where
    [(); WIDTH * HEIGHT]:,
{
    pub grid: [T; WIDTH * HEIGHT],
}

impl<T, const WIDTH: usize, const HEIGHT: usize> ModularArray<T, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
    T: Copy,
{
    pub fn new() -> Self
    where
        T: Default + Clone,
    {
        Self {
            grid: [T::default(); WIDTH * HEIGHT],
        }
    }
}

impl<T, const WIDTH: usize, const HEIGHT: usize> Index<(isize, isize)>
    for ModularArray<T, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
{
    type Output = T;

    fn index(&self, index: (isize, isize)) -> &Self::Output {
        let (x, y) = index;
        let x = x.rem_euclid(WIDTH as isize);
        let y = y.rem_euclid(HEIGHT as isize);
        self.grid.index(y as usize * WIDTH + x as usize)
    }
}
impl<T, const WIDTH: usize, const HEIGHT: usize> IndexMut<(isize, isize)>
    for ModularArray<T, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
{
    fn index_mut(&mut self, index: (isize, isize)) -> &mut Self::Output {
        let (x, y) = index;
        let x = x.rem_euclid(WIDTH as isize);
        let y = y.rem_euclid(HEIGHT as isize);
        self.grid.index_mut(y as usize * WIDTH + x as usize)
    }
}

impl<T, const WIDTH: usize, const HEIGHT: usize> Index<(usize, usize)>
    for ModularArray<T, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
{
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (x, y) = index;
        self.grid.index(y * WIDTH + x)
    }
}

impl<T, const WIDTH: usize, const HEIGHT: usize> IndexMut<(usize, usize)>
    for ModularArray<T, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
{
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (x, y) = index;
        self.grid.index_mut(y * WIDTH + x)
    }
}

impl<T: Default + Copy, const WIDTH: usize, const HEIGHT: usize> Default
    for ModularArray<T, WIDTH, HEIGHT>
where
    [(); WIDTH * HEIGHT]:,
{
    fn default() -> Self {
        Self::new()
    }
}
