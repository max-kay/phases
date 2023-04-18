use std::{
    fs::File,
    io::{BufWriter, Result, Write},
    sync::mpsc::{SendError, Sender},
    thread::JoinHandle,
};

#[derive(Clone)]
pub struct CsvLogger {
    sender: Sender<String>,
    columns: usize,
}

impl CsvLogger {
    pub fn new(
        path: String,
        header: String,
        categories: Vec<String>,
    ) -> (Self, JoinHandle<Result<()>>) {
        let (tx, rx) = std::sync::mpsc::channel::<String>();
        let columns = categories.len();
        let handle = std::thread::spawn(move || -> Result<()> {
            let file = File::create(path)?;
            let mut writer = BufWriter::new(file);
            writeln!(writer, "{}", header)?;
            writeln!(writer, "{}", categories.join(","))?;
            while let Ok(line) = rx.recv() {
                writeln!(writer, "{}", line)?;
            }
            writer.flush()
        });
        (
            Self {
                sender: tx,
                columns,
            },
            handle,
        )
    }

    pub fn send_row(&self, values: Vec<f32>) -> std::result::Result<(), SendError<String>> {
        assert_eq!(values.len(), self.columns);
        self.sender.send(
            values
                .iter()
                .map(|val| format!("{}", val))
                .collect::<Vec<String>>()
                .join(","),
        )
    }
}
