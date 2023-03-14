#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use gif::ExtensionData;
use num_integer::binomial;
use phases::{ArrayLatice, ModularArray, TriAtoms};
use rand_distr::Distribution;
use std::{fs::File, path::Path};

fn make_plot(ys: Vec<f32>, path: impl AsRef<Path>) {
    use plotters::prelude::*;

    let root = BitMapBackend::new(&path, (1280, 960)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    let max = *ys.iter().max_by(|y_1, y_2| y_1.total_cmp(y_2)).unwrap();
    let min = *ys.iter().min_by(|y_1, y_2| y_1.total_cmp(y_2)).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .caption("Energies", ("sans-serif", 30))
        .set_label_area_size(LabelAreaPosition::Left, 50)
        .set_label_area_size(LabelAreaPosition::Bottom, 50)
        .build_cartesian_2d(0.0..(ys.len() as f32), min..max)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(LineSeries::new(
            (0..ys.len()).zip(ys.iter()).map(|(x, y)| (x as f32, *y)),
            &BLACK,
        ))
        .unwrap();
}

#[derive(Clone, Copy)]
struct Distr<const N: u32>;

impl<const N: u32> Distribution<isize> for Distr<N> {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> isize {
        let val = rng.gen_range(0..(2_u32.pow(2_u32 * N)));
        let mut current_tot = 0;
        for i in 0..2 * N {
            current_tot += binomial(2 * N, i + 1);
            if val <= current_tot {
                return (i + 1) as isize - N as isize;
            }
        }
        unreachable!()
    }
}

fn energies(atom_1: TriAtoms, atom_2: TriAtoms) -> f32 {
    match (atom_1, atom_2) {
        (TriAtoms::A, TriAtoms::A) => 1.5,
        (TriAtoms::B, TriAtoms::B) => 0.5,
        (TriAtoms::C, TriAtoms::C) => 1.0,
        (TriAtoms::A, TriAtoms::B) | (TriAtoms::B, TriAtoms::A) => 0.3,
        (TriAtoms::A, TriAtoms::C) | (TriAtoms::C, TriAtoms::A) => -0.8,
        (TriAtoms::B, TriAtoms::C) | (TriAtoms::C, TriAtoms::B) => 1.0,
    }
}

const WIDTH: usize = 200;
const HEIGHT: usize = 200;
const STEPS: usize = 2_000_000;
const FRAMES: usize = 1000;
const LENGTH: usize = 1000; // in units of 10 ms

fn main() {
    let mut grid = ModularArray::<TriAtoms, WIDTH, HEIGHT>::new();

    for x in 50_usize..100 {
        for y in 60..110 {
            grid[(x, y)] = TriAtoms::B
        }
    }
    for x in 100_usize..WIDTH {
        for y in 60..110 {
            grid[(x, y)] = TriAtoms::C
        }
    }

    let mut latice =
        ArrayLatice::<_, WIDTH, HEIGHT>::new_from_grid(grid, energies, Some("my_seed"));

    let file = File::create("./out/swap_4.gif").expect("Error while creating file!");
    let palette: &[u8] = &[255, 127, 0, 0, 127, 255, 0, 255, 0];
    let mut encoder = gif::Encoder::new(file, WIDTH as u16, HEIGHT as u16, palette)
        .expect("Error while creating gif encoder");
    encoder
        .set_repeat(gif::Repeat::Infinite)
        .expect("Error while setting repeats!");
    encoder
        .write_extension(ExtensionData::new_control_ext(
            (LENGTH / FRAMES) as u16,
            gif::DisposalMethod::Any,
            true,
            None,
        ))
        .expect("Error while writing ExtensionData!");

    let start = std::time::Instant::now();
    let mut energies = Vec::with_capacity(STEPS);

    for i in 0..STEPS {
        latice.swap_dist_distr(Distr::<5>);
        energies.push(latice.tot_energy());
        if i % (STEPS / FRAMES) == 0 {
            let frame: &[u8] = (&latice.grid).into();
            let frame = gif::Frame::from_indexed_pixels(WIDTH as u16, HEIGHT as u16, frame, None);
            encoder
                .write_frame(&frame)
                .expect("Error while writing frame!");
        }
    }

    make_plot(energies, "out/swap_4.png");

    println!("took {:?}", std::time::Instant::now() - start);
}
