use crate::utils::{get_windows, make_bam_reader, Window};
use rust_htslib::bam::{Header, HeaderView, IndexedReader, Read, Record};

// we use these values to mark when the bam reader hasn't loaded the first reference/window.
const UNINIT_U32: u32 = u32::MAX - 1;
const UNINIT_USIZE: usize = usize::MAX - 1;
const UNINIT_I64: i64 = i64::MAX - 1;

pub struct WindowedBamReader {
    reader: IndexedReader,
    window_size: Option<i64>,
    pub raw_header: Header,
    meta_header: HeaderView,
    pub windows: Vec<Vec<Window>>,
    pub cur_window: Window,
    cur_window_idx: usize,
    pub cur_ref: u32,
}

impl WindowedBamReader {
    pub fn new(file_name: &String, num_threads: usize, window_size: Option<i64>) -> Self {
        let (raw_header, reader) = make_bam_reader(file_name, num_threads);
        let meta_header = reader.header().clone();
        let cur_ref = UNINIT_U32;

        Self {
            reader,
            window_size,
            raw_header,
            meta_header,
            windows: vec![],
            cur_window: Window {
                start: UNINIT_I64,
                end: UNINIT_I64,
            },
            cur_window_idx: usize::MAX,
            cur_ref,
        }
    }

    pub fn window_records(&mut self) -> impl Iterator<Item = Record> + '_ {
        self.reader.records().flatten().take_while(|record| {
            record.pos() >= self.cur_window.start && record.pos() < self.cur_window.end
        })
    }

    pub fn next_reference(&mut self) -> bool {
        self.cur_ref = match self.cur_ref {
            UNINIT_U32 => 0,
            _ => self.cur_ref + 1,
        };

        self.reader.fetch((self.cur_ref, 0, u32::MAX)).unwrap();

        if self.cur_ref >= self.meta_header.target_count() {
            false
        } else {
            let windows = get_windows(
                self.window_size,
                self.meta_header.target_len(self.cur_ref).unwrap() as i64,
            )
            .chunks(3)
            .map(|chunk| chunk.to_vec())
            .collect::<Vec<Vec<Window>>>();

            self.windows = windows;
            self.cur_window_idx = usize::MAX;
            true
        }
    }

    pub fn next_window(&mut self) -> bool {
        self.cur_window_idx = match self.cur_window_idx {
            UNINIT_USIZE => 0,
            _ => self.cur_window_idx + 1,
        };

        if let Some(window_chunk) = self.windows.get(self.cur_window_idx) {
            self.cur_window.start = window_chunk.first().unwrap().start;
            self.cur_window.end = window_chunk.last().unwrap().end;
            true
        } else {
            false
        }
    }
}
