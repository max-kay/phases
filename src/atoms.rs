use crate::MyRng;
use rand::Rng;
use std::{hash::Hash, ops::Deref, usize};

pub trait RandAtom: Default + Eq + PartialEq + Hash + Deref<Target = u8> {
    type Concentration: Copy;
    fn vacancy() -> Self;
    fn with_concentration(rng: &mut MyRng, cs: Self::Concentration) -> Self;
}
pub trait Mark: RandAtom {
    /// # Safety
    /// this function is allowed to manipulate bits of the underlying data
    /// this has to be undone by using the `.unmark()` function
    unsafe fn mark(&mut self);
    fn unmark(&mut self);
    fn is_marked(&self) -> bool;
}

#[derive(Hash, PartialEq, Eq, Default, Clone, Copy)]
#[repr(transparent)]
pub struct BinAtom(u8);
#[derive(Clone, Copy)]
pub struct BinConcentration(f64);

impl BinConcentration {
    pub fn new(c_a: f64, c_b: f64) -> Self {
        Self(c_a / (c_a + c_b))
    }

    pub fn get_c_a(&self) -> f64 {
        self.0
    }
}
impl BinAtom {
    pub fn new(num: u8) -> Self {
        assert!(num <= 1);
        Self(num)
    }
}

impl Deref for BinAtom {
    type Target = u8;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl RandAtom for BinAtom {
    type Concentration = BinConcentration;

    fn vacancy() -> Self {
        Self(0b0100)
    }

    fn with_concentration(rng: &mut MyRng, cs: Self::Concentration) -> Self {
        if rng.gen_bool(cs.0) {
            Self(0b0000)
        } else {
            Self(0b0001)
        }
    }
}

impl Mark for BinAtom {
    unsafe fn mark(&mut self) {
        self.0 |= 0b1000_0000
    }

    fn unmark(&mut self) {
        self.0 &= 0b0111_1111
    }

    fn is_marked(&self) -> bool {
        (self.0 & 0b1000_0000) != 0
    }
}

pub trait Energies<A: RandAtom> {
    fn get_interaction_energy(&self, a_1: A, a_2: A) -> f32;
    fn as_dict(&self) -> String;
}

impl Energies<BinAtom> for [f32; 4] {
    fn get_interaction_energy(&self, a_1: BinAtom, a_2: BinAtom) -> f32 {
        // Safety this is save as with the last & the index is ensured to be < 4
        // this is also correct with the vacancy value of 0b0000_0100 but it depends on the
        // concrete implementation of the move vacancy method
        unsafe { *self.get_unchecked((((*a_1 << 1) + *a_2) & 0b0000_0011) as usize) }
    }

    fn as_dict(&self) -> String {
        format!(
            "{{(0, 0): {}, (0, 1): {}, (1, 1): {}}}",
            self[0], self[1], self[3]
        )
    }
}
