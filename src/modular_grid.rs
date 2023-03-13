use std::ops::{Index, IndexMut};

#[derive(Clone)]
pub struct ModularGrid<T> {
    width: usize,
    height: usize,
    pub grid: Vec<T>,
}

impl<T> ModularGrid<T> {
    pub fn new(width: usize, height: usize) -> Self
    where
        T: Default + Clone,
    {
        Self {
            width,
            height,
            grid: vec![T::default(); width * height],
        }
    }

    pub fn modularize(width: usize, height: usize, grid: Vec<T>) -> Self {
        assert!(
            grid.len() == width * height,
            "len of grid was {}, but should be width*height: {}",
            grid.len(),
            width * height
        );
        Self {
            width,
            height,
            grid,
        }
    }
}

impl<T> Index<(isize, isize)> for ModularGrid<T> {
    type Output = T;

    fn index(&self, index: (isize, isize)) -> &Self::Output {
        let (x, y) = index;
        let x = x.rem_euclid(self.width as isize);
        let y = y.rem_euclid(self.height as isize);
        self.grid.index(y as usize * self.width + x as usize)
    }
}
impl<T> IndexMut<(isize, isize)> for ModularGrid<T> {
    fn index_mut(&mut self, index: (isize, isize)) -> &mut Self::Output {
        let (x, y) = index;
        let x = x.rem_euclid(self.width as isize);
        let y = y.rem_euclid(self.height as isize);
        self.grid.index_mut(y as usize * self.width + x as usize)
    }
}

impl<T> Index<(usize, usize)> for ModularGrid<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (x, y) = index;
        self.grid.index(y * self.width + x)
    }
}

impl<T> IndexMut<(usize, usize)> for ModularGrid<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (x, y) = index;
        self.grid.index_mut(y * self.width + x)
    }
}


