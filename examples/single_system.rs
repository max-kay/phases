#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use gif::Frame;
use phases::{
    anim::prepare_encoder,
    plots::{float_plot},
    Lattice,
};

// model parameters
const N_ATOMS: usize = 2;
type Atom = phases::Atom<N_ATOMS>;
type Concentration = phases::Concentration<N_ATOMS>;
const WIDTH: usize = 400;
const HEIGHT: usize = 400;
const STEPS: usize = WIDTH * HEIGHT * 400;

fn energies(a1: Atom, a2: Atom) -> f32 {
    match (*a1, *a2) {
        (0, 0) => -1.0,
        (1, 1) => -1.0,
        (0, 1) | (1, 0) => 3.0,
        // (2, 2) => 0.0,
        // (2, 0) | (0, 2) => 0.0,
        // (2, 1) | (1, 2) => 0.0,
        _ => panic!(),
    }
}

// temperature
const START: f32 = 0.002;
const RISE: f32 = 0.000002;
fn beta(i: usize) -> f32 {
    START * (RISE * i as f32).exp()
    // 30.0
}

// gif
const FRAMES: usize = 100;
const LENGTH: usize = 5000; // in ms

fn main() {
    let name = "mmm";

    let mut lattice = Lattice::<N_ATOMS, WIDTH, HEIGHT>::new(
        energies,
        Some("my_seed"),
        Some(Concentration::new([1.0, 1.0])),
    );

    let mut encoder = prepare_encoder(
        format!("./out/{}.gif", name),
        WIDTH as u16,
        HEIGHT as u16,
        Some((LENGTH / FRAMES) as u16),
    );

    let start = std::time::Instant::now();
    let mut energies = Vec::with_capacity(STEPS);
    for i in 0..STEPS {
        lattice.monte_carlo_swap(beta(i));
        energies.push(lattice.internal_energy());
        if i % (STEPS / FRAMES) == 0 {
            encoder
                .write_frame(&Frame::from_indexed_pixels(
                    WIDTH as u16,
                    HEIGHT as u16,
                    lattice.as_bytes(),
                    None,
                ))
                .expect("Error while writing frame!");
        }
    }
    println!("took {:?}", std::time::Instant::now() - start);

    let mut encoder = prepare_encoder(
        format!("./out/{}_last_frame.gif", name),
        WIDTH as u16,
        HEIGHT as u16,
        None,
    );
    encoder
        .write_frame(&Frame::from_indexed_pixels(
            WIDTH as u16,
            HEIGHT as u16,
            lattice.as_bytes(),
            None,
        ))
        .expect("Error while writing frame!");

    let steps_per_latice: Vec<f32> = (0..STEPS).map(|i| i as f32 / WIDTH as f32 /HEIGHT as f32).collect();
    float_plot(
        steps_per_latice.clone(),
        energies.clone(),
        format!("out/{}.png", name),
        "",
        "steps",
        "internal energy",
    )
    .unwrap();
    float_plot(
        steps_per_latice,
        (0..STEPS).map(|i| 1.0 / beta(i)).collect(),
        format!("out/{}_temp_curve.png", name),
        "Temperature curve",
        "steps",
        "temperature",
    )
    .unwrap();
    float_plot(
        (0..STEPS).map(|i| 1.0 / beta(i)).collect(),
        energies,
        format!("out/{}_temp.png", name),
        "",
        "temperature",
        "internal energy",
    )
    .unwrap();
}
