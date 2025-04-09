pub type Ngram<'a> = &'a str;

// for repeatedly processing UMIs of the same length, we can store the chunk size to avoid
// calculating it each time.
#[derive(Debug)]
pub struct NgramMaker {
    chunk_size: usize,
    pub num_chunks: usize,
    string_len: usize,
}

impl<'a> NgramMaker {
    pub fn new(num_chunks: usize, string_len: usize) -> Self {
        let mut chunk_size = string_len / num_chunks;

        if string_len % num_chunks != 0 {
            chunk_size += 1;
        }

        Self {
            chunk_size,
            num_chunks,
            string_len,
        }
    }

    pub fn ngrams_to_vec(&mut self, string: &'a str, out_vec: &mut Vec<&'a str>) {
        // the output vector needs to have placeholders
        // TODO: this is still slightly slower than the old method for some reason, find out why.
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
    let mut ngram_maker = NgramMaker::new(2, s.len());
    let mut ngrams = vec!["NILL"; ngram_maker.num_chunks];
    // 15 / 2 = 7
    // 7 + 1 = 8
    ngram_maker.ngrams_to_vec(s, &mut ngrams);
    println!("{:?}", &ngrams);
    assert!(ngrams == vec!["AGCTCTAG".to_string(), "CTACGAG".to_string()])
}

#[test]
fn test_split2() {
    let s = "GGGGGGGCCCCCCGGGGGGCCCCCC";
    let mut ngram_maker = NgramMaker::new(3, s.len());
    let mut ngrams = vec!["NILL"; ngram_maker.num_chunks];
    // 25 / 3 = 8
    // 8 + 1 = 9
    ngram_maker.ngrams_to_vec(s, &mut ngrams);
    println!("{:?}", ngrams);
    assert!(
        ngrams
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
    let mut ngram_maker = NgramMaker::new(2, s.len());
    let mut ngrams = vec!["NULL"; ngram_maker.num_chunks];
    // 6 / 2 = 3
    ngram_maker.ngrams_to_vec(s, &mut ngrams);
    println!("{:?}", ngrams);
    assert!(*ngrams == vec!["GTC".to_string(), "TAC".to_string()])
}
