use std::{fs::File, io::Write};

use chrono::Utc;
use phases::{
    anim::prepare_encoder, get_energies_dict, logs::CsvLogger, run_python, Array2d, Stats, System,
};

// model parameters
type Atom = phases::NumAtom<2>;
type Concentration = phases::NumC<2>;
const WIDTH: usize = 512;
const HEIGHT: usize = 512;
const STEPS: usize = WIDTH * HEIGHT * 1000;

// temperature
const START: f32 = 50.0;
const END: f32 = 0.01;
fn temp(i: usize) -> f32 {
    START * ((END / START).ln() / STEPS as f32 * i as f32).exp()
}

// gif
const FRAMES: usize = 60;
const LENGTH: usize = 2000; // in ms

// logs
const LOG_ENTRIES: usize = 1000;

fn main() {
    let start = std::time::Instant::now();

    let concentration = Concentration::new([1.0, 1.0]);

    let name = format!("b_2_t_{}", Utc::now().format("%Y-%m-%d_%H-%M"));

    make_system_file(&name, concentration).unwrap();

    let path = format!("out/logs/{}.csv", name);

    let mut categories = vec!["step".to_owned(), "temp".to_owned(), "energy".to_owned()];
    categories.append(&mut Stats::get_categories(Some("atom 0 ")));
    categories.append(&mut Stats::get_categories(Some("atom 2 ")));

    let (logger, handle) = CsvLogger::new(
        path,
        "file generated as log to github maxkay/phases".to_owned(),
        categories,
    );

    let mut encoder = prepare_encoder(
        format!("out/gifs/{}.gif", name),
        WIDTH as u16,
        HEIGHT as u16,
        Some((LENGTH / FRAMES) as u16),
    );

    let mut system =
        System::<Array2d<Atom, WIDTH, HEIGHT>>::new(energies, Some("my_seed"), Some(concentration));

    for i in 0..STEPS {
        system.move_vacancy(1.0 / temp(i));
        if i % (STEPS / LOG_ENTRIES) == 0 {
            let mut values = vec![
                i as f32 / (WIDTH * HEIGHT) as f32,
                temp(i),
                system.internal_energy() / (WIDTH * HEIGHT) as f32,
            ];
            let stats = system.get_region_stats();
            values.append(&mut stats[&Atom::new(0)].as_f32_vec());
            values.append(&mut stats[&Atom::new(1)].as_f32_vec());
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

    let mut encoder = prepare_encoder(
        format!("out/gifs/{}_last.gif", name),
        WIDTH as u16,
        HEIGHT as u16,
        None,
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
    let energies_dict = get_energies_dict(energies);
    writeln!(file, "{name}")?;
    writeln!(file, "Energies")?;
    writeln!(file, "{energies_dict}")?;
    writeln!(file, "Width, Height")?;
    writeln!(file, "{WIDTH}, {HEIGHT}")?;
    writeln!(file, "Steps per Lattice Site")?;
    writeln!(file, "{}", STEPS as f32 / (WIDTH * HEIGHT) as f32)?;
    writeln!(file, "Concentration")?;
    writeln!(file, "{:?}", concentration.get_cs())?;
    Ok(())
}

#[inline(always)]
fn energies(a1: Atom, a2: Atom) -> f32 {
    // Safety: this is safe because NumAtom<2> can only be 0 or 1
    // and thus the shift is equal to multiplying by 2
    unsafe { *[-4.0, 3.0, 3.0, -1.0].get_unchecked(((*a1 << 1) + *a2) as usize) }
}
