use phases::{Array2d, System};

// model parameters
const N_ATOMS: usize = 2;
type Atom = phases::NumAtom<N_ATOMS>;
type Concentration = phases::NumC<N_ATOMS>;
const WIDTH: usize = 200;
const HEIGHT: usize = 200;
const STEPS: usize = WIDTH * HEIGHT * 500;

// temperature
const START: f32 = 50.0;
const END: f32 = 0.01;
fn temp(i: usize) -> f32 {
    START * ((END / START).ln() / STEPS as f32 * i as f32).exp()
}

fn main() {
    let start = std::time::Instant::now();
    let concentration = Concentration::new([1.0, 1.0]);
    let mut system =
        System::<Array2d<Atom, WIDTH, HEIGHT>>::new(energies, Some("my_seed"), Some(concentration));
    for i in 0..STEPS {
        system.move_vacancy(1.0 / temp(i));
    }
    println!("finshed running model, took: {:?}", start.elapsed());
}

#[inline(always)]
fn energies(a1: Atom, a2: Atom) -> f32 {
    [-4.0, 3.0, 3.0, -1.0][(*a1 * 2 + *a2) as usize]
}
