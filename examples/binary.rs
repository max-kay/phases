#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use std::{fs::File, io::Write, sync::atomic::AtomicU64};

use chrono::Utc;
use phases::{Array2d, System};
use rayon::prelude::*;

// model parameters
type Atom = phases::NumAtom<2>;
type Concentration = phases::NumC<2>;
const WIDTH: usize = 512;
const HEIGHT: usize = 512;
const STEPS: usize = WIDTH * HEIGHT * 300;
const EQUILIBRIUM_STEPS: usize = WIDTH * HEIGHT * 100;

fn energies(a1: Atom, a2: Atom) -> f32 {
    match (*a1, *a2) {
        (0, 0) => -4.0,
        (1, 1) => -1.0,
        (0, 1) | (1, 0) => 3.0,
        _ => panic!(),
    }
}

// temp
const TEMP_STEPS: u32 = 100;
const START_TEMP: f32 = 300.0;

// concentration
const CONCENTRATION_STEPS: usize = 28;

static PROGRESS_COUNTER: AtomicU64 = AtomicU64::new(1);

fn main() {
    let start = std::time::Instant::now();

    let temps: Vec<f32> = (0..TEMP_STEPS)
        .map(|i| START_TEMP / TEMP_STEPS as f32 * i as f32)
        .rev()
        .collect();
    let concentrations: Vec<f64> = (0..CONCENTRATION_STEPS)
        .map(|i| (i as f64 / CONCENTRATION_STEPS as f64) * 0.8 + 0.1)
        .collect();

    let results: Vec<(Vec<f32>, Vec<f32>)> = concentrations
        .par_iter()
        .map(|c_a| {
            run_model_with_concentration(Concentration::new([*c_a, 1.0 - c_a]), temps.clone())
        })
        .collect();

    let file_name = format!("logs_bin/data_{}.csv", Utc::now().format("%Y-%m-%d_%H-%M"));
    let mut file = File::create(&file_name).expect("error while creating file");
    writeln!(file, "c,temp,energy,heat capacity").unwrap();

    for (c, (int_energies, heat_capacities)) in concentrations.iter().zip(results.iter()) {
        for ((t, int_energy), heat_capacity) in temps
            .iter()
            .zip(int_energies.iter())
            .zip(heat_capacities.iter())
        {
            writeln!(file, "{:?},{:?},{:?},{:?}", c, t, int_energy, heat_capacity).unwrap();
        }
    }

    println!("took {:?}", start.elapsed());
}

fn run_model_with_concentration(
    concentration: Concentration,
    temps: Vec<f32>,
) -> (Vec<f32>, Vec<f32>) {
    // let mut encoder = prepare_encoder(
    //     format!("out/{:.0}.gif", concentration.get_cs()[0] * 100.0),
    //     WIDTH as u16,
    //     HEIGHT as u16,
    //     Some((5000 / TEMP_STEPS) as u16),
    // );

    let mut avg_int_energies: Vec<f32> = Vec::new();
    let mut heat_capacity: Vec<f32> = Vec::new();
    let mut lattice =
        System::<Array2d<Atom, WIDTH, HEIGHT>>::new(energies, None, Some(concentration));
    for temp in temps {
        let beta = 1.0 / temp;

        for _ in 0..EQUILIBRIUM_STEPS {
            lattice.move_vacancy(beta);
        }

        let mut int_energies: Vec<f32> = Vec::with_capacity(STEPS);
        for _ in 0..STEPS {
            lattice.move_vacancy(beta);
            int_energies.push(lattice.internal_energy())
        }

        // doing it like this because of numerics
        let avg_energy: f32 = int_energies.iter().sum::<f32>() / int_energies.len() as f32;
        let variance = int_energies
            .iter()
            .map(|e| (*e - avg_energy).powi(2))
            .sum::<f32>()
            / int_energies.len() as f32;
        let avg_energy = avg_energy / (WIDTH * HEIGHT) as f32;
        let variance = variance / (WIDTH * HEIGHT) as f32 / (WIDTH * HEIGHT) as f32;
        avg_int_energies.push(avg_energy);
        heat_capacity.push(variance / temp.powi(2));
        let progress = PROGRESS_COUNTER.fetch_add(1, std::sync::atomic::Ordering::AcqRel);
        println!(
            "{} of {}",
            progress,
            TEMP_STEPS * CONCENTRATION_STEPS as u32
        );
        // encoder
        //     .write_frame(&Frame::from_indexed_pixels(
        //         WIDTH as u16,
        //         HEIGHT as u16,
        //         lattice.as_bytes(),
        //         None,
        //     ))
        //     .unwrap();
    }
    (avg_int_energies, heat_capacity)
}
