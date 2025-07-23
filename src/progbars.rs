use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

pub struct ProgressTracker {
    _prog: Option<MultiProgress>,
    ref_bar: Option<ProgressBar>,
    window_bar: Option<ProgressBar>,
    pub coord_bar: ProgressBar,
}

impl ProgressTracker {
    pub fn initialize_main(num_references: u32, show: bool) -> Self {
        if !show {
            Self {
                _prog: None,
                ref_bar: None,
                window_bar: None,
                coord_bar: ProgressBar::hidden(),
            }
        } else {
            let _prog = MultiProgress::new();
            _prog.clear().ok();
            let ref_bar = _prog.add(make_reference_bar(num_references as u64));
            ref_bar.set_prefix("REFERENCE");

            let ref_bar = Some(ref_bar);

            let window_bar = _prog.add(make_windowbar(0)); // Length set in initialize_windows
            window_bar.set_prefix("WINDOW");
            window_bar.tick();

            let window_bar = Some(window_bar);

            let coord_bar = _prog.add(make_coordbar(1));
            coord_bar.set_prefix("COORDINATES");
            coord_bar.tick();

            let coord_bar = coord_bar;

            let _prog = Some(_prog);

            Self {
                _prog,
                ref_bar,
                window_bar,
                coord_bar,
            }
        }
    }

    pub fn initialize_windows(&mut self, num_windows: usize) {
        if let Some(window_bar) = &self.window_bar {
            window_bar.reset();
            window_bar.set_length(num_windows as u64);
            window_bar.tick();
        }
    }

    pub fn intake_reads_msg(&mut self) {
        if let Some(window_bar) = &self.window_bar {
            window_bar.set_message("Intaking reads...");
        }
    }

    pub fn update_window_reads(&mut self, num_reads: u64) {
        if let Some(window_bar) = &self.window_bar {
            window_bar.set_message(format!("{} reads in window", num_reads));
            window_bar.tick();
        }
    }

    pub fn next_ref(&mut self) {
        if let (Some(ref_bar), Some(window_bar)) = (&self.ref_bar, &self.window_bar) {
            ref_bar.inc(1);
            window_bar.reset();
            ref_bar.tick();
        }
    }

    pub fn next_window(&mut self) {
        if let Some(window_bar) = &self.window_bar {
            window_bar.inc(1);
            window_bar.tick();
            self.coord_bar.reset();
        }
    }

    pub fn finish(&self) {
        if let (Some(window_bar), Some(ref_bar)) = (&self.window_bar, &self.ref_bar) {
            ref_bar.finish_with_message("DONE");
            window_bar.finish_with_message("DONE");
            self.coord_bar.finish_and_clear();
        }
    }
}

pub fn make_reference_bar(num_refs: u64) -> ProgressBar {
    let pb = ProgressBar::new(num_refs);
    pb.set_style(
        ProgressStyle::with_template(
            "{prefix:<15} {human_pos:>3}/{human_len:<3} {msg:<15} {spinner}",
        )
        .unwrap()
        .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(100)); // Auto-tick
    pb
}

pub fn make_windowbar(num_windows: u64) -> ProgressBar {
    let pb = ProgressBar::new(num_windows);
    pb.set_style(
        ProgressStyle::with_template(
            "{prefix:<15} {human_pos:>3}/{human_len:<3} {msg:<15} {spinner}",
        )
        .unwrap()
        .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(100)); // Auto-tick
    pb
}

pub fn make_coordbar(num_coords: u64) -> ProgressBar {
    let pb = ProgressBar::new(num_coords);
    pb.set_style(
        ProgressStyle::with_template(
            "{prefix:<15} {human_pos:>3}/{human_len:<3} {bar:40.cyan/blue}",
        )
        .unwrap(),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(100));
    pb
}
