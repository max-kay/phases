#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use phases::{ArrayLatice, TerAtoms, TerConcentration};
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

fn energies(atom_1: TerAtoms, atom_2: TerAtoms) -> f32 {
    match (atom_1, atom_2) {
        (TerAtoms::A, TerAtoms::A) => -2.0,
        (TerAtoms::B, TerAtoms::B) => -1.8,
        (TerAtoms::C, TerAtoms::C) => 0.1,
        (TerAtoms::A, TerAtoms::B) | (TerAtoms::B, TerAtoms::A) => 5.0,
        (TerAtoms::A, TerAtoms::C) | (TerAtoms::C, TerAtoms::A) => -1.0,
        (TerAtoms::B, TerAtoms::C) | (TerAtoms::C, TerAtoms::B) => -1.0,
    }
}

const WIDTH: usize = 200;
const HEIGHT: usize = 200;
const STEPS: usize = 200_000_000;
const FRAMES: usize = 100;

fn main() {
    let name = "mmdm";

    let mut latice = ArrayLatice::<_, WIDTH, HEIGHT>::new(
        energies,
        Some("my_seed"),
        Some(TerConcentration::new(0.6, 0.3, 0.1)),
    );

    let file = File::create(format!("./out/{}.gif", name)).expect("Error while creating file!");
    let palette: &[u8] = &[255, 127, 0, 0, 127, 255, 0, 255, 0];
    let mut encoder = gif::Encoder::new(file, WIDTH as u16, HEIGHT as u16, palette)
        .expect("Error while creating gif encoder");
    encoder
        .set_repeat(gif::Repeat::Infinite)
        .expect("Error while setting repeats!");
    // encoder
    //     .write_extension(ExtensionData::new_control_ext(
    //         (LENGTH / FRAMES) as u16,
    //         gif::DisposalMethod::Any,
    //         true,
    //         None,
    //     ))
    //     .expect("Error while writing ExtensionData!");

    let start = std::time::Instant::now();
    let mut energies = Vec::with_capacity(STEPS);
    let mut temperature = 0.2;
    for i in 0..STEPS {
        latice.monte_carlo_swap(temperature);
        energies.push(latice.tot_energy());
        if i % (STEPS / FRAMES) == 0 {
            temperature *= 0.99_f32.powi(10);
            let frame = gif::Frame::from_indexed_pixels(
                WIDTH as u16,
                HEIGHT as u16,
                latice.as_bytes(),
                None,
            );
            encoder
                .write_frame(&frame)
                .expect("Error while writing frame!");
        }
    }

    let file =
        File::create(format!("./out/{}_last_frame.gif", name)).expect("Error while creating file!");
    let mut encoder = gif::Encoder::new(file, WIDTH as u16, HEIGHT as u16, palette)
        .expect("Error while creating gif encoder");
    encoder
        .set_repeat(gif::Repeat::Infinite)
        .expect("Error while setting repeats!");

    let frame = gif::Frame::from_indexed_pixels(WIDTH as u16, HEIGHT as u16, latice.as_bytes(), None);
    encoder
        .write_frame(&frame)
        .expect("Error while writing frame!");

    make_plot(energies, format!("out/{}.png", name));

    println!("took {:?}", std::time::Instant::now() - start);
}
