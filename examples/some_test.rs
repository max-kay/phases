use std::{fs::File, io::Write, sync::atomic::AtomicU64};

use chrono::Utc;
use phases::{
    logs::CsvLogger, run_python, BinAtom as Atom, BinConcentration as Concentration, Energies,
    FastArray, StreamingStats, System,
};
use rayon::prelude::*;

// model parameters
const SIDE: usize = 128;
const POW: usize = 7;
const STEPS: usize = SIDE * SIDE * 4000;
const FIRST_STEPS: usize = SIDE * SIDE * 4000;
const EQUILIBRIUM_STEPS: usize = SIDE * SIDE * 10_000;

const ENERGIES: [f32; 4] = [-1.0, -0.75, -0.75, -1.0];

// temp
const TEMP_STEPS: usize = 15;
const START_TEMP: f32 = 2.0;

// concentration
const CONCENTRATION_STEPS: usize = 8;

static PROGRESS_COUNTER: AtomicU64 = AtomicU64::new(1);

fn main() {
    let start = std::time::Instant::now();
    let name = format!("some_test_{}", Utc::now().format("%Y-%m-%d_%H-%M"));

    make_system_file(&name).unwrap();

    let temps: Vec<f32> = (0..TEMP_STEPS)
        .map(|i| START_TEMP / TEMP_STEPS as f32 * i as f32)
        .rev()
        .collect();
    let concentrations: Vec<f64> = (0..CONCENTRATION_STEPS)
        .map(|i| (i as f64 / (CONCENTRATION_STEPS - 1) as f64) * 0.92 + 0.04)
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
                Concentration::new(*c_a, 1.0 - c_a),
                temps.clone(),
                logger.clone(),
            )
        })
        .collect();
    handle
        .join()
        .expect("failed to log")
        .expect("failed to log");
    println!("finished running {}, took {:?}", name, start.elapsed());

    run_python("python/b_f.py", &name)
}

fn run_model_with_concentration(concentration: Concentration, temps: Vec<f32>, logger: CsvLogger) {
    let mut system = System::<FastArray<Atom, SIDE, POW>, _>::new(ENERGIES, None, concentration);
    for _ in 0..FIRST_STEPS {
        system.move_vacancy(1.0 / temps[0]);
    }
    for temp in temps {
        let beta = 1.0 / temp;

        for _ in 0..EQUILIBRIUM_STEPS {
            system.move_vacancy(beta);
        }

        let mut stats = StreamingStats::new();
        for _ in 0..STEPS {
            system.move_vacancy(beta);
            stats.add_value(system.internal_energy())
        }

        logger
            .send_row(vec![
                concentration.get_c_a() as f32,
                temp,
                stats.avg() / (SIDE * SIDE) as f32,
                stats.variance() / (temp * temp) / (SIDE * SIDE) as f32,
            ])
            .unwrap();

        let progress = PROGRESS_COUNTER.fetch_add(1, std::sync::atomic::Ordering::AcqRel);
        println!("{} of {}", progress, TEMP_STEPS * CONCENTRATION_STEPS);
    }
    std::mem::drop(logger)
}

fn make_system_file(name: &String) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(format!("out/systems/{}.txt", name))?;
    let energies_dict = ENERGIES.as_dict();

    writeln!(file, "{}", name)?;
    writeln!(file, "energies:")?;
    writeln!(file, "{}", energies_dict)?;
    writeln!(file, "width, height")?;
    writeln!(file, "{}, {}", SIDE, SIDE)?;
    writeln!(file, "steps_per_site: for equilibrium, for measurement")?;
    writeln!(
        file,
        "{},{}",
        EQUILIBRIUM_STEPS / SIDE / SIDE,
        STEPS / SIDE / SIDE
    )?;
    writeln!(file, "start temp, temp steps")?;
    writeln!(file, "{}, {}", START_TEMP, TEMP_STEPS)?;
    writeln!(file, "concentration steps")?;
    writeln!(file, "{}", CONCENTRATION_STEPS)?;
    Ok(())
}
