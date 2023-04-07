use std::{fs::File, path::Path};

use gif::{ExtensionData, Repeat};

const PALETTE: &[u8] = &[
    0xE6, 0x9F, 0x00, 0x56, 0xB4, 0xE9, 0x00, 0x9E, 0x73, 0xF0, 0xE4, 0x42, 0x00, 0x72, 0xB2, 0xD5,
    0x5E, 0x00, 0xCC, 0x79, 0xA7, 0x00, 0x00, 0x00,
];

pub fn prepare_encoder(
    path: impl AsRef<Path>,
    width: u16,
    height: u16,
    mus_per_frame: Option<u16>,
) -> gif::Encoder<File> {
    let file = File::create(path).expect("Error while creating file!");
    let mut encoder =
        gif::Encoder::new(file, width, height, PALETTE).expect("Error while creating gif encoder");
    encoder
        .set_repeat(Repeat::Infinite)
        .expect("Error while setting repeats!");
    if let Some(mus_per_frame) = mus_per_frame {
        encoder
            .write_extension(ExtensionData::new_control_ext(
                mus_per_frame,
                gif::DisposalMethod::Any,
                true,
                None,
            ))
            .expect("Error while writing ExtensionData!");
    }
    encoder
}
