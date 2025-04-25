use smol_str::SmolStr;
use std::cell::{RefCell, RefMut};

#[derive(Debug, Default)]
/// A utility struct to generate a fixed number of ngrams for a given string. It maintains an internal, fixed-length
/// vector for holding generated ngrams, avoiding excessive allocations.
///
/// Ngram generation is fastest when dealing with an alphabet of same-size words; in case a word is
/// longer than the initial size used to set up the [NgramMaker], the last ngram will hold all remaining characters.
///
/// To use it to generate ngrams for connecting strings >= K edits apart, it should be initialized with
/// num_chunks = K + 1.
pub struct NgramMaker {
    chunk_size: usize,
    out_vec: RefCell<Vec<SmolStr>>,
    num_chunks: usize,
}

impl NgramMaker {
    pub fn new(num_chunks: usize, string_len: usize) -> Self {
        let mut chunk_size = string_len / num_chunks;

        if string_len % num_chunks != 0 {
            chunk_size += 1;
        }

        let out_vec = RefCell::new(vec![SmolStr::new("NILL"); num_chunks]);

        Self {
            chunk_size,
            out_vec,
            num_chunks,
        }
    }

    pub fn ngrams(&self, s: &str) -> RefMut<Vec<SmolStr>> {
        self.ngrams_to_ref(s, self.out_vec.borrow_mut());
        self.out_vec.borrow_mut()
    }

    fn ngrams_to_ref(&self, string: &str, mut out_vec: RefMut<Vec<SmolStr>>) {
        let mut start = 0;
        let mut end = 0;
        let s_len = string.len();

        let mut cur_idx = 0;

        while start < s_len && cur_idx < self.num_chunks {
            end = (start + self.chunk_size).min(s_len);
            out_vec[cur_idx] = SmolStr::new(&string[start..end]);
            start = end;
            cur_idx += 1;
        }

        // we subtract 1 from end because it's non-inclusive
        let rem = s_len.saturating_sub(1).saturating_sub(end - 1);

        // if we have remaining characters, add them to last chunk
        // this is used for strings that don't fit into n [Self::num_chunks] when divided by
        // [Self::chunk_size]
        if rem > 0 {
            let last = out_vec[self.num_chunks - 1].as_str();
            let last = SmolStr::new([last, &string[end..]].concat());
            out_vec[self.num_chunks - 1] = last;
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

#[test]
fn test_split4() {
    let s1 = "GTCTAC";
    let s2 = "GTCTACG";

    let ngram_maker = NgramMaker::new(2, s1.len());
    // we init our ngram maker with enough num_chunks to split s1 in two pieces, but s2 is longer
    // by 1 char.
    // we expect to have the second ngram of s2 to hold the remaining char.
    // 6 / 2 = 3
    let ngrams = ngram_maker.ngrams(s1);
    println!("{:?}", ngrams);
    assert!(*ngrams == vec!["GTC".to_string(), "TAC".to_string()]);

    drop(ngrams);

    let ngrams = ngram_maker.ngrams(s2);
    println! {"{:?}", ngrams};
    assert!(*ngrams == vec!["GTC".to_string(), "TACG".to_string()]);
}
