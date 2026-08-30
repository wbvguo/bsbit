use bsbit_hts::TextRecordLimits;

pub const fn text_limits() -> TextRecordLimits {
    TextRecordLimits::new(512, 8, 128, 256, 384, 1_024, 384)
}

pub fn buffer_capacity(control: u8) -> usize {
    usize::from(control % 31 + 1)
}
