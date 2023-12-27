use std::{fs::File, io::Write};

use chrono::Utc;
use phases::{
    anim::{self, prepare_file_encoder},
    logs::CsvLogger,
    run_python, BinAtom as Atom, BinConcentration as Concentration, ClusterStats, Energies,
    FastArray, System,
};
// model parameters
const SIDE: usize = 256;
const POW: usize = 8;
const STEPS: usize = SIDE * SIDE * 40_000;
const ENERGIES: [f32; 4] = [-1.0, -0.75, -0.75, -1.0];

// temperature
const START: f32 = 8.0;
const END: f32 = 0.01;
fn temp(i: usize) -> f32 {
    START * ((END / START).ln() / STEPS as f32 * i as f32).exp()
}

// gif
const FRAMES: usize = 60;
const LENGTH: usize = 2000; // in ms

// logs
const LOG_ENTRIES: usize = 10_000;

fn main() {
    let start = std::time::Instant::now();

    let concentration = Concentration::new(1.0, 1.0);

    let name = format!("b_2_t_{}", Utc::now().format("%Y-%m-%d_%H-%M"));

    make_system_file(&name, concentration).unwrap();

    let path = format!("out/logs/{}.csv", name);

    let mut categories = vec!["step".to_owned(), "temp".to_owned(), "energy".to_owned()];
    categories.append(&mut ClusterStats::get_categories(Some("atom 0 ")));
    categories.append(&mut ClusterStats::get_categories(Some("atom 2 ")));

    let (logger, handle) = CsvLogger::new(
        path,
        "file generated as log to github maxkay/phases".to_owned(),
        categories,
    );

    let mut encoder =
        prepare_file_encoder(
            format!("out/gifs/{}.gif", name),
            SIDE as u16,
            SIDE as u16,
            Some((LENGTH / FRAMES) as u16),
            anim::PALETTE,
        );

    let mut system =
        System::<FastArray<Atom, SIDE, POW>, _>::new(ENERGIES, Some("my_seed"), concentration);

    for i in 0..STEPS {
        system.move_vacancy(1.0 / temp(i));
        if i % (STEPS / LOG_ENTRIES) == 0 {
            let mut values = vec![
                i as f32 / (SIDE * SIDE) as f32,
                temp(i),
                system.internal_energy() / (SIDE * SIDE) as f32,
            ];
            for distr in system.count_all_clusters() {
                values.append(&mut ClusterStats::from_map_atom(distr).as_vec_f32());
            }
            logger.send_row(values).expect("error while sending row");
        }
        if i % (STEPS / FRAMES) == 0 {
            let frame = system.get_frame();
            encoder
                .write_frame(&frame)
                .expect("Error while writing frame!");
        }
    }

    std::mem::drop(logger);
    if let Err(err) = handle.join() {
        eprintln!("an error occurred while logging: {:?}", err)
    };

    let mut encoder = prepare_file_encoder(
        format!("out/gifs/{}_last.gif", name),
        SIDE as u16,
        SIDE as u16,
        None,
        anim::PALETTE,
    );
    encoder
        .write_frame(&system.get_frame())
        .expect("Error while writing frame!");

    println!("took {:?}", start.elapsed());

    run_python("python/b_t.py", &name)
}

fn make_system_file(
    name: &String,
    concentration: Concentration,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(format!("out/systems/{}.txt", name))?;
    let energies_dict = ENERGIES.as_dict();
    writeln!(file, "{name}")?;
    writeln!(file, "Energies")?;
    writeln!(file, "{energies_dict}")?;
    writeln!(file, "Width, Height")?;
    writeln!(file, "{SIDE}, {SIDE}")?;
    writeln!(file, "Steps per Lattice Site")?;
    writeln!(file, "{}", STEPS as f32 / (SIDE * SIDE) as f32)?;
    writeln!(file, "Concentration")?;
    writeln!(file, "{:?}", concentration.get_c_a())?;
    Ok(())
}
