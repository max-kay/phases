use std::{fs::File, io::Write};

use chrono::Utc;
use gif::Frame;
use phases::{
    anim::prepare_encoder, get_energies_dict, logs::CsvLogger, run_python, Array2d, System,
};

// model parameters
const N_ATOMS: usize = 2;
type Atom = phases::NumAtom<N_ATOMS>;
type Concentration = phases::NumC<N_ATOMS>;
const WIDTH: usize = 150;
const HEIGHT: usize = 200;
const STEPS: usize = WIDTH * HEIGHT * 500;

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
    let categories = vec!["temp".to_owned(), "energy".to_owned()];

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
            logger
                .send_row(vec![temp(i), system.internal_energy()])
                .expect("error while sending row");
        }
        if i % (STEPS / FRAMES) == 0 {
            encoder
                .write_frame(&Frame::from_indexed_pixels(
                    WIDTH as u16,
                    HEIGHT as u16,
                    system.as_bytes(),
                    None,
                ))
                .expect("Error while writing frame!");
        }
    }
    println!("finshed running model, took: {:?}", start.elapsed());

    let mut encoder = prepare_encoder(
        format!("out/gifs/{}_last.gif", name),
        WIDTH as u16,
        HEIGHT as u16,
        None,
    );
    encoder
        .write_frame(&Frame::from_indexed_pixels(
            WIDTH as u16,
            HEIGHT as u16,
            system.as_bytes(),
            None,
        ))
        .expect("Error while writing frame!");

    std::mem::drop(logger);
    if let Err(err) = handle.join() {
        eprintln!("an error occurred while logging: {:?}", err)
    };
    println!("took {:?}", start.elapsed());

    run_python("python/b_t.py", &name)
}

fn make_system_file(
    name: &String,
    concentration: phases::NumC<N_ATOMS>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(format!("out/systems/{}.txt", name))?;
    let energies_dict = get_energies_dict(energies);

    writeln!(file, "energies")?;
    writeln!(file, "{}", energies_dict)?;
    writeln!(file, "width, height")?;
    writeln!(file, "{},{}", WIDTH, HEIGHT)?;
    writeln!(file, "steps")?;
    writeln!(file, "{}", STEPS)?;
    writeln!(file, "{:?}", concentration.get_cs())?;
    Ok(())
}

#[inline(always)]
fn energies(a1: Atom, a2: Atom) -> f32 {
    [-4.0, 3.0, 3.0, -1.0][(*a1 * 2 + *a2) as usize]
}
