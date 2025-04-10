pub type Ngram<'a> = &'a str;
use std::cell::{Ref, RefCell, RefMut};

#[derive(Debug, Default)]
/// A utility struct to generate ngrams for a given string of a given, fixed length. It maintains an internal, fixed-length
/// vector for holding generated ngrams, avoiding excessive allocations.
///
/// To use it to generate ngrams for connecting strings >= K edits apart, it should be initialized with
/// num_chunks = K + 1.
pub struct NgramMaker<'a> {
    chunk_size: usize,
    pub num_chunks: usize,
    string_len: usize,
    out_vec: RefCell<Vec<Ngram<'a>>>,
}

impl<'a> NgramMaker<'a> {
    pub fn new(num_chunks: usize, string_len: usize) -> Self {
        let mut chunk_size = string_len / num_chunks;

        if string_len % num_chunks != 0 {
            chunk_size += 1;
        }

        let out_vec = RefCell::new(vec!["NILL"; num_chunks]);

        Self {
            chunk_size,
            num_chunks,
            string_len,
            out_vec,
        }
    }

    pub fn ngrams(&self, s: &'a str) -> Ref<Vec<&'a str>> {
        self.ngrams_to_ref(s, self.out_vec.borrow_mut());
        self.out_vec.borrow()
    }

    fn ngrams_to_ref(&self, string: &'a str, mut out_vec: RefMut<Vec<&'a str>>) {
        let mut start = 0;

        let mut cur_idx = 0;
        while start < self.string_len {
            let end = (start + self.chunk_size).min(self.string_len);
            out_vec[cur_idx] = &string[start..end];
            start = end;
            cur_idx += 1;
        }
    }
}

#[test]
fn test_split1() {
    let s = "AGCTCTAGCTACGAG";
    let ngram_maker = NgramMaker::new(2, s.len());
    // 15 / 2 = 7
    // 7 + 1 = 8
    let ngrams = ngram_maker.ngrams(s);
    println!("{:?}", &ngrams);
    assert!(*ngrams == vec!["AGCTCTAG", "CTACGAG"])
}

#[test]
fn test_split2() {
    let s = "GGGGGGGCCCCCCGGGGGGCCCCCC";
    let ngram_maker = NgramMaker::new(3, s.len());
    // 25 / 3 = 8
    // 8 + 1 = 9
    let ngrams = ngram_maker.ngrams(s);
    println!("{:?}", ngrams);
    assert!(
        *ngrams
            == vec![
                "GGGGGGGCC".to_string(),
                "CCCCGGGGG".to_string(),
                "GCCCCCC".to_string()
            ]
    )
}

#[test]
fn test_split3() {
    let s = "GTCTAC";
    let ngram_maker = NgramMaker::new(2, s.len());
    // 6 / 2 = 3
    let ngrams = ngram_maker.ngrams(s);
    println!("{:?}", ngrams);
    assert!(*ngrams == vec!["GTC".to_string(), "TAC".to_string()])
}
