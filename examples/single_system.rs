#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
use gif::Frame;
use phases::{anim::prepare_encoder, plots::float_plot, Array2d, System};
use rand_distr::Distribution;

// model parameters
const N_ATOMS: usize = 3;
type Atom = phases::NumAtom<N_ATOMS>;
type Concentration = phases::NumC<N_ATOMS>;
const WIDTH: usize = 400;
const HEIGHT: usize = 400;
const STEPS: usize = WIDTH * HEIGHT * 4000;

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

// temperature
const START: f32 = 50.0;
const END: f32 = 0.1;
fn temp(i: usize) -> f32 {
    START * ((END / START).ln() / STEPS as f32 * i as f32).exp()
    // 50.0
}

// gif
const FRAMES: usize = 100;
const LENGTH: usize = 5000; // in ms

fn main() {
    let name = "meug";

    let mut lattice = System::<Array2d<Atom, WIDTH, HEIGHT>>::new(
        energies,
        Some("my_seed"),
        Some(Concentration::new([1.0, 0.3, 1.0])),
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
        lattice.monte_carlo_swap_distr(MyDistr, 1.0 / temp(i));
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

    let steps_per_latice: Vec<f32> = (0..STEPS)
        .map(|i| i as f32 / WIDTH as f32 / HEIGHT as f32)
        .collect();
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
        (0..STEPS).map(temp).collect(),
        format!("out/{}_temp_curve.png", name),
        "Temperature curve",
        "steps",
        "temperature",
    )
    .unwrap();
    float_plot(
        (0..STEPS).map(temp).collect(),
        energies,
        format!("out/{}_temp.png", name),
        "",
        "temperature",
        "internal energy",
    )
    .unwrap();
    println!("took {:?}", std::time::Instant::now() - start);
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
