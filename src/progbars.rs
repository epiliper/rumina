use indicatif::{ProgressBar, ProgressStyle};

pub fn make_readbar() -> ProgressBar {
    ProgressBar::new_spinner().with_style(ProgressStyle::with_template("{msg}").unwrap())
}

pub fn make_windowbar(num_windows: u64) -> ProgressBar {
    ProgressBar::with_style(
        ProgressBar::new(num_windows),
        ProgressStyle::with_template(
            "{prefix}:\t\t       {human_pos}/{human_len:7} {msg:15}{spinner.white}",
        )
        .unwrap(),
    )
}

pub fn make_coordbar(num_coords: u64) -> ProgressBar {
    ProgressBar::with_style(
        ProgressBar::new(num_coords),
        ProgressStyle::with_template("{prefix}:\t{human_pos}/{human_len:7} {bar:40.cyan/blue}")
            .unwrap(),
    )
}
