use plotters::prelude::*;
use std::{error::Error, path::Path};

const WIDTH: u32 = 2048;
const HEIGHT: u32 = 1536;

pub fn index_plot(
    ys: Vec<f32>,
    path: impl AsRef<Path>,
    caption: &str,
    x_desc: &str,
    y_desc: &str,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new(&path, (WIDTH, HEIGHT)).into_drawing_area();
    root.fill(&WHITE)?;

    let max = *ys
        .iter()
        .max_by(|y_1, y_2| y_1.total_cmp(y_2))
        .expect("the data should be non empty");
    let min = *ys
        .iter()
        .min_by(|y_1, y_2| y_1.total_cmp(y_2))
        .expect("the data should be non empty");

    let mut chart = ChartBuilder::on(&root)
        .margin(HEIGHT / 15)
        .caption(caption, ("sans-serif", HEIGHT / 20))
        .set_label_area_size(LabelAreaPosition::Left, HEIGHT / 30)
        .set_label_area_size(LabelAreaPosition::Bottom, HEIGHT / 30)
        .build_cartesian_2d(0.0..(ys.len() as f32), min..max)?;

    chart
        .configure_mesh()
        .disable_mesh()
        .bold_line_style(BLACK)
        .x_label_style(("sans_serif", HEIGHT / 40))
        .x_label_formatter(&|x| format!("{:.0}", x))
        .x_desc(x_desc)
        .y_label_style(("sans_serif", HEIGHT / 40))
        .y_label_formatter(&|y| format!("{}", y))
        .y_desc(y_desc)
        .draw()?;

    chart.draw_series(
        (0..ys.len())
            .zip(ys.iter())
            .map(|(x, y)| Circle::new((x as f32, *y), 0.005, BLACK.filled())),
    )?;
    Ok(())
}

pub fn float_plot(
    xs: Vec<f32>,
    ys: Vec<f32>,
    path: impl AsRef<Path>,
    caption: &str,
    x_desc: &str,
    y_desc: &str,
) -> Result<(), Box<dyn Error>> {
    assert_eq!(xs.len(), ys.len(), "xs and ys must be of same length");

    let root = BitMapBackend::new(&path, (WIDTH, HEIGHT)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_x = *xs
        .iter()
        .max_by(|x_1, x_2| x_1.total_cmp(x_2))
        .expect("the data should be non empty");
    let min_x = *xs
        .iter()
        .min_by(|x_1, x_2| x_1.total_cmp(x_2))
        .expect("the data should be non empty");

    let max_y = *ys
        .iter()
        .max_by(|y_1, y_2| y_1.total_cmp(y_2))
        .expect("the data should be non empty");
    let min_y = *ys
        .iter()
        .min_by(|y_1, y_2| y_1.total_cmp(y_2))
        .expect("the data should be non empty");

    let mut chart = ChartBuilder::on(&root)
        .margin(HEIGHT / 15)
        .caption(caption, ("sans-serif", HEIGHT / 20))
        .set_label_area_size(LabelAreaPosition::Left, HEIGHT / 30)
        .set_label_area_size(LabelAreaPosition::Bottom, HEIGHT / 30)
        .build_cartesian_2d(min_x..max_x, min_y..max_y)?;

    chart
        .configure_mesh()
        .disable_mesh()
        .bold_line_style(BLACK)
        .x_label_style(("sans_serif", HEIGHT / 40))
        .x_label_formatter(&|x| format!("{}", x))
        .x_desc(x_desc)
        .y_label_style(("sans_serif", HEIGHT / 40))
        .y_label_formatter(&|y| format!("{}", y))
        .y_desc(y_desc)
        .draw()?;

    chart.draw_series(
        xs.iter()
            .zip(ys.iter())
            .map(|(x, y)| Circle::new((*x, *y), 0.005, BLACK.filled())),
    )?;
    Ok(())
}
pub fn line_plot(
    xs: Vec<f32>,
    ys: Vec<f32>,
    path: impl AsRef<Path>,
    caption: &str,
    x_desc: &str,
    y_desc: &str,
) -> Result<(), Box<dyn Error>> {
    assert_eq!(xs.len(), ys.len(), "xs and ys must be of same length");

    let root = BitMapBackend::new(&path, (WIDTH, HEIGHT)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_x = *xs
        .iter()
        .max_by(|x_1, x_2| x_1.total_cmp(x_2))
        .expect("the data should be non empty");
    let min_x = *xs
        .iter()
        .min_by(|x_1, x_2| x_1.total_cmp(x_2))
        .expect("the data should be non empty");

    let max_y = *ys
        .iter()
        .max_by(|y_1, y_2| y_1.total_cmp(y_2))
        .expect("the data should be non empty");
    let min_y = *ys
        .iter()
        .min_by(|y_1, y_2| y_1.total_cmp(y_2))
        .expect("the data should be non empty");

    let mut chart = ChartBuilder::on(&root)
        .margin(HEIGHT / 15)
        .caption(caption, ("sans-serif", HEIGHT / 20))
        .set_label_area_size(LabelAreaPosition::Left, HEIGHT / 30)
        .set_label_area_size(LabelAreaPosition::Bottom, HEIGHT / 30)
        .build_cartesian_2d(min_x..max_x, min_y..max_y)?;

    chart
        .configure_mesh()
        .disable_mesh()
        .bold_line_style(BLACK)
        .x_label_style(("sans_serif", HEIGHT / 40))
        .x_label_formatter(&|x| format!("{}", x))
        .x_desc(x_desc)
        .y_label_style(("sans_serif", HEIGHT / 40))
        .y_label_formatter(&|y| format!("{}", y))
        .y_desc(y_desc)
        .draw()?;

    chart.draw_series(LineSeries::new(
        xs.iter().zip(ys.iter()).map(|(x, y)| (*x, *y)),
        BLACK,
    ))?;
    Ok(())
}
