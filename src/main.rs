#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use phases::{ArrayLatice, TriAtoms};
use std::fs::File;

fn energies(atom_1: TriAtoms, atom_2: TriAtoms) -> f32 {
    match (atom_1, atom_2) {
        (TriAtoms::A, TriAtoms::A) => -0.2,
        (TriAtoms::B, TriAtoms::B) => 0.5,
        (TriAtoms::C, TriAtoms::C) => 5.0,
        (TriAtoms::A, TriAtoms::B) | (TriAtoms::B, TriAtoms::A) => 0.3,
        (TriAtoms::A, TriAtoms::C) | (TriAtoms::C, TriAtoms::A) => -0.2,
        (TriAtoms::B, TriAtoms::C) | (TriAtoms::C, TriAtoms::B) => 1.0,
    }
}

const WIDTH: usize = 200;
const HEIGHT: usize = 200;
const STEPS: usize = 1000000;
const FRAMES: usize = 100;

fn main() {
    let mut grid =
        ArrayLatice::<_, WIDTH, HEIGHT>::new(energies, Some("my_seed"), Some((0.1, 0.5)));

    let file = File::create("./out/test.gif").expect("Error while creating file!");
    let palette: &[u8] = &[255, 127, 0, 0, 127, 255, 0, 255, 0];
    let mut encoder = gif::Encoder::new(file, WIDTH as u16, HEIGHT as u16, palette)
        .expect("Error while creating gif encoder");
    encoder
        .set_repeat(gif::Repeat::Finite(0))
        .expect("Error while setting repeats!");

    let start = std::time::Instant::now();

    for i in 0..STEPS {
        grid.swap_rand();
        if i % (STEPS / FRAMES) == 0 {
            let frame: &[u8] = (&grid.grid).into();
            let frame = gif::Frame::from_indexed_pixels(WIDTH as u16, HEIGHT as u16, frame, None);
            encoder
                .write_frame(&frame)
                .expect("Error while writing frame!");
        }
    }

    println!("took {:?}", std::time::Instant::now() - start);
}
