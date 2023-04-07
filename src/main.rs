#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use std::io::Write;
use std::{fs::File, sync::atomic::AtomicU64};

use chrono::Utc;
use phases::Lattice;
use rayon::prelude::*;

// model parameters
type Atom = phases::Atom<2>;
type Concentration = phases::Concentration<2>;
const WIDTH: usize = 128;
const HEIGHT: usize = 128;
const STEPS: usize = WIDTH * HEIGHT * 100;
const EQUILIBRIUM_STEPS: usize = WIDTH * HEIGHT * 200;

fn energies(a1: Atom, a2: Atom) -> f32 {
    match (*a1, *a2) {
        (0, 0) => -2.0,
        (1, 1) => -1.0,
        (0, 1) | (1, 0) => 3.0,
        // (2, 2) => 0.0,
        // (2, 0) | (0, 2) => 0.0,
        // (2, 1) | (1, 2) => 0.0,
        _ => panic!(),
    }
}

// temp
const TEMP_STEPS: u32 = 30;
const START_TEMP: f32 = 50.0;

// concentration
const CONCENTRATION_STEPS: usize = 21;

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
        .map(|c_a| run_model_at_concentration(Concentration::new([*c_a, 1.0 - c_a]), temps.clone()))
        .collect();
    let mut file = File::create(format!(
        "logs/data_{}.csv",
        Utc::now().format("%Y-%m-%d_%H-%M")
    ))
    .expect("error while creating file");
    writeln!(
        file,
        "concentration a,temperature,internal energy U,heat capacity"
    )
    .unwrap();
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

fn run_model_at_concentration(
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
    let mut lattice = Lattice::<2, WIDTH, HEIGHT>::new(energies, Some("seed"), Some(concentration));
    for temp in temps {
        let beta = 1.0 / temp;

        for _ in 0..EQUILIBRIUM_STEPS {
            lattice.monte_carlo_swap(beta);
        }

        let mut cumulative_int_energy = 0.0;
        let mut cumulative_int_energy_squared = 0.0;

        for _ in 0..STEPS {
            lattice.monte_carlo_swap(beta);
            let int_energy = lattice.internal_energy();
            cumulative_int_energy += int_energy;
            cumulative_int_energy_squared += int_energy * int_energy;
        }

        let avg_energy = cumulative_int_energy / STEPS as f32 / (WIDTH as f32 * HEIGHT as f32);
        let avg_energy_sq = cumulative_int_energy_squared
            / STEPS as f32
            / (WIDTH as f32 * WIDTH as f32 * HEIGHT as f32 * HEIGHT as f32);
        avg_int_energies.push(avg_energy);
        heat_capacity.push(dbg!(avg_energy_sq - avg_energy * avg_energy) / temp.powi(2));
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
