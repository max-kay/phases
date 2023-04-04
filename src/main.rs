#![allow(incomplete_features, dead_code, unused_imports)]
#![feature(generic_const_exprs)]
use gif::{ExtensionData, Frame, Repeat};
use phases::ArrayLatice;
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

// model parameters
const N_ATOMS: usize = 2;
type Atom = phases::Atom<N_ATOMS>;
type Concentration = phases::Concrete<N_ATOMS>;
const WIDTH: usize = 500;
const HEIGHT: usize = 500;
const STEPS: usize = WIDTH * HEIGHT * 1000;

fn energies(a1: Atom, a2: Atom) -> f32 {
    match (*a1, *a2) {
        (0, 0) => 0.0,
        (0, 1) | (1, 0) => 0.0,
        (1, 1) => 0.0,
        _ => panic!(),
    }
}

// temperature
const START: f32 = 0.002;
const RISE: f32 = 0.0000002;
fn beta(i: usize) -> f32 {
    START * (RISE * i as f32).exp()
}

const FRAMES: usize = 100;
const LENGTH: usize = 5000; // in ms

fn main() {
    println!("start setup");
    let start = std::time::Instant::now();

    let name = "mmm";

    // binary
    let mut latice = ArrayLatice::<N_ATOMS, WIDTH, HEIGHT>::new(energies, Some("my_seed"), None);
    let palette: &[u8] = &[25, 127, 0, 0, 200, 180];

    // // ternary
    // let palette: &[u8] = &[255, 127, 0, 25, 127, 255, 0, 255, 60];

    let file = File::create(format!("./out/{}.gif", name)).expect("Error while creating file!");
    let mut encoder = gif::Encoder::new(file, WIDTH as u16, HEIGHT as u16, palette)
        .expect("Error while creating gif encoder");
    encoder
        .set_repeat(Repeat::Infinite)
        .expect("Error while setting repeats!");
    encoder
        .write_extension(ExtensionData::new_control_ext(
            (LENGTH / FRAMES) as u16,
            gif::DisposalMethod::Any,
            true,
            None,
        ))
        .expect("Error while writing ExtensionData!");
    println!("took {:?}", std::time::Instant::now() - start);

    let start = std::time::Instant::now();
    let mut energies = Vec::with_capacity(STEPS);
    for i in 0..STEPS {
        latice.monte_carlo_swap(beta(i));
        energies.push(latice.internal_energy());
        if i % (STEPS / FRAMES) == 0 {
            // println!("{} of {}", i / (STEPS / FRAMES) + 1, FRAMES);
            encoder
                .write_frame(&Frame::from_indexed_pixels(
                    WIDTH as u16,
                    HEIGHT as u16,
                    latice.as_bytes(),
                    None,
                ))
                .expect("Error while writing frame!");
        }
    }

    let file =
        File::create(format!("./out/{}_last_frame.gif", name)).expect("Error while creating file!");
    let mut encoder = gif::Encoder::new(file, WIDTH as u16, HEIGHT as u16, palette)
        .expect("Error while creating gif encoder");
    encoder
        .set_repeat(Repeat::Infinite)
        .expect("Error while setting repeats!");

    println!("took {:?}", std::time::Instant::now() - start);

    encoder
        .write_frame(&Frame::from_indexed_pixels(
            WIDTH as u16,
            HEIGHT as u16,
            latice.as_bytes(),
            None,
        ))
        .expect("Error while writing frame!");

    make_plot(energies, format!("out/{}.png", name));
}
