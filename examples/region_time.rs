use std::fs::File;

use phases::{logs::CsvLogger, Array2d, NumAtom, RegionCounter, Stats};

const PATH: &str = "out/gifs/b_2_t_2023-04-18_15-41.gif";
const SIZE: usize = 512;

fn main() {
    let (logger, handle) = CsvLogger::new(
        "out/test.csv".to_owned(),
        "".to_owned(),
        vec![
            "min".to_owned(),
            "quart_1".to_owned(),
            "median".to_owned(),
            "quart_3".to_owned(),
            "max".to_owned(),
        ],
    );
    let file = File::open(PATH).unwrap();
    let decoder = gif::DecodeOptions::new();
    let mut decoder = decoder.read_info(file).unwrap();
    let mut tot_time = std::time::Duration::from_secs(0);
    while let Some(frame) = decoder.read_next_frame().unwrap() {
        let mut arr = Box::new([[0_u8; SIZE]; SIZE]);
        for i in 0..SIZE {
            for j in 0..SIZE {
                arr[i][j] = frame.buffer[i * SIZE + j];
            }
        }
        let mut grid: Array2d<NumAtom<2>, SIZE, SIZE> = unsafe { std::mem::transmute(arr) };
        let start = std::time::Instant::now();
        let stats = Stats::gen_stats(grid.count_regions());
        let stats = &stats[&NumAtom::<2>::new(0)];
        logger.send_row(stats.as_f32_vec()).unwrap();

        tot_time += start.elapsed();
    }
    std::mem::drop(logger);
    handle.join().unwrap().unwrap();
    println!("took: {:?} to count the regions", tot_time)
}
