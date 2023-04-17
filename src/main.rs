#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use std::{
    fs::File,
    io::Write,
    sync::{atomic::AtomicU64, mpsc::Sender},
};

use chrono::Utc;
use phases::{Array2d, System};
use rand_distr::Distribution;
use rayon::prelude::*;

// model parameters
type Atom = phases::NumAtom<3>;
type Concentration = phases::NumC<3>;
const WIDTH: usize = 256;
const HEIGHT: usize = 512;
const STEPS: usize = WIDTH * HEIGHT * 200;
const EQUILIBRIUM_STEPS: usize = WIDTH * HEIGHT * 300;

fn energies(a1: Atom, a2: Atom) -> f32 {
    match (*a1, *a2) {
        (0, 0) => -1.0,
        (1, 1) => -1.0,
        (0, 1) | (1, 0) => 3.0,
        (2, 2) => 0.0,
        (2, 0) | (0, 2) => -2.0,
        (2, 1) | (1, 2) => 0.5,
        _ => unreachable!(),
    }
}

// temp
const TEMP_STEPS: u32 = 30;
const START_TEMP: f32 = 100.0;

// concentration
const CONCENTRATION_STEPS: usize = 21;
const BLUE_CONCENTRATIONS: &[f64] = &[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];

static PROGRESS_COUNTER: AtomicU64 = AtomicU64::new(1);

fn main() {
    let start = std::time::Instant::now();
    let (tx, rx) = std::sync::mpsc::channel::<String>();
    std::thread::spawn(move || {
        let file_name = format!("out/logs/data_{}.csv", Utc::now().format("%Y-%m-%d_%H-%M"));
        let mut file = File::create(&file_name).expect("error while creating file");
        writeln!(file,
            "A simulation of three different atom types with bonding energies 00: {}, 01: {}, 02: {}, 11: {}, 12: {}, 22: {}",
            energies(Atom::new(0), Atom::new(0)),
            energies(Atom::new(0), Atom::new(1)),
            energies(Atom::new(0), Atom::new(2)),
            energies(Atom::new(1), Atom::new(1)),
            energies(Atom::new(1), Atom::new(2)),
            energies(Atom::new(2), Atom::new(2)),
        ).expect("error while writing header");

        writeln!(file, "c0,c1,c2,temp,energy,heat capacity").unwrap();

        while let Ok(line) = rx.recv() {
            writeln!(file, "{}", line).expect("an error occured while wrinting an line");
        }
    });

    for c_blue in BLUE_CONCENTRATIONS {
        let temps: Vec<f32> = (0..TEMP_STEPS)
            .map(|i| START_TEMP / TEMP_STEPS as f32 * i as f32)
            .rev()
            .collect();
        let a_concentrations: Vec<f64> = (0..CONCENTRATION_STEPS)
            .map(|i| (i as f64 / CONCENTRATION_STEPS as f64) * 0.8 + 0.1)
            .collect();

        let _: Vec<_> = a_concentrations
            .par_iter()
            .map_with(tx.clone(), |tx, c_a| {
                let concentration = Concentration::new([
                    (1.0 - c_blue) * c_a,
                    *c_blue,
                    (1.0 - c_blue) * (1.0 - c_a),
                ]);
                run_model_at_concentration(concentration, temps.clone(), tx.clone())
            })
            .collect();
    }
    println!("took {:?}", start.elapsed());
}

fn run_model_at_concentration(
    concentration: Concentration,
    temps: Vec<f32>,
    sender: Sender<String>,
) {
    let mut lattice =
        System::<Array2d<Atom, WIDTH, HEIGHT>>::new(energies, None, Some(concentration));
    for temp in temps {
        let beta = 1.0 / temp;

        for _ in 0..EQUILIBRIUM_STEPS {
            lattice.monte_carlo_swap_distr(MyDistr, beta);
        }

        let mut int_energies: Vec<f32> = Vec::with_capacity(STEPS);
        for _ in 0..STEPS {
            lattice.monte_carlo_swap_distr(MyDistr, beta);
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
        let variance = variance / (WIDTH * HEIGHT * WIDTH * HEIGHT) as f32;

        let progress = PROGRESS_COUNTER.fetch_add(1, std::sync::atomic::Ordering::AcqRel);

        sender
            .send(format!(
                "{:?},{:?},{:?},{:?},{:?},{:?}",
                concentration.get_cs()[0],
                concentration.get_cs()[1],
                concentration.get_cs()[2],
                temp,
                avg_energy,
                variance / temp / temp,
            ))
            .expect("an error occured while sending string");
        println!(
            "{} of {}",
            progress,
            TEMP_STEPS * CONCENTRATION_STEPS as u32 * BLUE_CONCENTRATIONS.len() as u32
        );
    }
}

#[derive(Copy, Clone)]
pub struct MyDistr;

impl Distribution<(isize, isize)> for MyDistr {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> (isize, isize) {
        (
            if rng.gen() { 1 } else { -1 },
            if rng.gen() { 1 } else { -1 },
        )
    }
}
