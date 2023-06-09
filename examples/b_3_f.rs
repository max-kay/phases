use std::{fs::File, io::Write, mem::drop, sync::atomic::AtomicU64};

use chrono::Utc;
use phases::{logs::CsvLogger, run_python, Array3d, Energies, System};
use rayon::prelude::*;

// model parameters
type Atom = phases::NumAtom<2>;
type Concentration = phases::NumC<2>;
const WIDTH: usize = 64;
const HEIGHT: usize = 64;
const DEPTH: usize = 64;
const STEPS: usize = WIDTH * HEIGHT * DEPTH * 100;
const EQUILIBRIUM_STEPS: usize = WIDTH * HEIGHT * DEPTH * 60;

const ENERGIES: [f32; 4] = [-1.0, -0.75, -0.75, -1.0];

// temp
const TEMP_STEPS: usize = 100;
const START_TEMP: f32 = 150.0;

// concentration
const CONCENTRATION_STEPS: usize = 30;

static PROGRESS_COUNTER: AtomicU64 = AtomicU64::new(1);

fn main() {
    let start = std::time::Instant::now();
    let name = format!("b_3_f_{}", Utc::now().format("%Y-%m-%d_%H-%M"));

    make_system_file(&name).unwrap();

    let temps: Vec<f32> = (0..TEMP_STEPS)
        .map(|i| START_TEMP / TEMP_STEPS as f32 * i as f32)
        .rev()
        .collect();
    let concentrations: Vec<f64> = (0..CONCENTRATION_STEPS)
        .map(|i| (i as f64 / CONCENTRATION_STEPS as f64) * 0.8 + 0.1)
        .collect();

    let (logger, handle) = CsvLogger::new(
        format!("out/logs/{}.csv", name),
        "file generated as log to github maxkay/phases".to_owned(),
        vec![
            "c".to_owned(),
            "temp".to_owned(),
            "energy".to_owned(),
            "heat capacity".to_owned(),
        ],
    );

    let _: Vec<_> = concentrations
        .par_iter()
        .map_with(logger, |logger, c_a| {
            run_model_with_concentration(
                Concentration::new([*c_a, 1.0 - c_a]),
                temps.clone(),
                logger.clone(),
            )
        })
        .collect();
    handle
        .join()
        .expect("failed to log")
        .expect("failed to log");
    println!("took {:?}", start.elapsed());

    run_python("python/b_f.py", &name)
}

fn run_model_with_concentration(concentration: Concentration, temps: Vec<f32>, logger: CsvLogger) {
    let mut system =
        System::<Array3d<Atom, WIDTH, HEIGHT, DEPTH>, _>::new(ENERGIES, None, concentration);
    for temp in temps {
        let beta = 1.0 / temp;

        for _ in 0..EQUILIBRIUM_STEPS {
            system.move_vacancy(beta);
        }

        let mut int_energies: Vec<f32> = Vec::with_capacity(STEPS);
        for _ in 0..STEPS {
            system.move_vacancy(beta);
            int_energies.push(system.internal_energy())
        }

        // doing it like this because of numerics
        let avg_energy: f32 = int_energies.iter().sum::<f32>() / int_energies.len() as f32;
        let variance = int_energies
            .iter()
            .map(|e| (*e - avg_energy).powi(2))
            .sum::<f32>()
            / int_energies.len() as f32;
        let avg_energy = avg_energy / (WIDTH * HEIGHT * DEPTH) as f32;
        let variance = variance / (WIDTH * HEIGHT * DEPTH) as f32 / (WIDTH * HEIGHT * DEPTH) as f32;

        logger
            .send_row(vec![
                concentration.get_cs()[0] as f32,
                temp,
                avg_energy,
                variance / temp / temp,
            ])
            .unwrap();

        let progress = PROGRESS_COUNTER.fetch_add(1, std::sync::atomic::Ordering::AcqRel);
        println!("{} of {}", progress, TEMP_STEPS * CONCENTRATION_STEPS);
    }
    drop(logger)
}

fn make_system_file(name: &String) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(format!("out/systems/{}.txt", name))?;
    let energies_dict = ENERGIES.as_dict();

    writeln!(file, "{}", name)?;
    writeln!(file, "energies:")?;
    writeln!(file, "{}", energies_dict)?;
    writeln!(file, "width, height, depth")?;
    writeln!(file, "{},{},{}", WIDTH, HEIGHT, DEPTH)?;
    writeln!(file, "steps_per_site: for equilibrium, for measurement")?;
    writeln!(
        file,
        "{},{}",
        EQUILIBRIUM_STEPS / WIDTH / HEIGHT / DEPTH,
        STEPS / WIDTH / HEIGHT / DEPTH
    )?;
    writeln!(file, "start temp, temp steps")?;
    writeln!(file, "{}, {}", START_TEMP, TEMP_STEPS)?;
    writeln!(file, "concentration steps")?;
    writeln!(file, "{}", CONCENTRATION_STEPS)?;
    Ok(())
}
