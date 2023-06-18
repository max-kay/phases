use std::{fs::File, path::Path};

use gif::{ExtensionData, Repeat};

pub const PALETTE: &[u8] = &[
    0xE6, 0x9F, 0x00, 0x56, 0xB4, 0xE9, 0x00, 0x9E, 0x73, 0xF0, 0xE4, 0x42, 0x00, 0x72, 0xB2, 0xD5,
    0x5E, 0x00, 0xCC, 0x79, 0xA7, 0x00, 0x00, 0x00,
];

pub fn prepare_file_encoder(
    path: impl AsRef<Path>,
    width: u16,
    height: u16,
    hundredths_per_frame: Option<u16>,
    palette: &[u8],
) -> gif::Encoder<File> {
    let file = File::create(path).expect("Error while creating file!");
    let mut encoder =
        gif::Encoder::new(file, width, height, palette).expect("Error while creating gif encoder");
    encoder
        .set_repeat(Repeat::Infinite)
        .expect("Error while setting repeats!");
    if let Some(hundredths_per_frame) = hundredths_per_frame {
        encoder
            .write_extension(ExtensionData::new_control_ext(
                hundredths_per_frame,
                gif::DisposalMethod::Any,
                true,
                None,
            ))
            .expect("Error while writing ExtensionData!");
    }
    encoder
}

pub fn prepare_vec_encoder(
    width: u16,
    height: u16,
    hundredths_per_frame: Option<u16>,
    palette: &[u8],
) -> gif::Encoder<Vec<u8>> {
    let mut encoder = gif::Encoder::new(Vec::new(), width, height, palette)
        .expect("Error while creating gif location");
    encoder
        .set_repeat(Repeat::Infinite)
        .expect("Error while setting repeats!");
    if let Some(hundredths_per_frame) = hundredths_per_frame {
        encoder
            .write_extension(ExtensionData::new_control_ext(
                hundredths_per_frame,
                gif::DisposalMethod::Any,
                true,
                None,
            ))
            .expect("Error while writing ExtensionData!");
    }
    encoder
}
