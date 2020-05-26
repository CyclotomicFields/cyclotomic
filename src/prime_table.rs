use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::Read;

type ZPlus = u64;

pub struct PrimeTableReader<'a> {
    directory_path: &'a Path
}

impl<'a> PrimeTableReader<'a> {
    pub fn new(directory_path: &Path) -> PrimeTableReader {
        assert!(directory_path.is_dir());
        return PrimeTableReader { directory_path };
    }

    pub fn first_million_primes(&self) -> Vec<ZPlus> {
        let primes_text_file_path: PathBuf = self.directory_path.join("primes1.txt");
        assert!(primes_text_file_path.is_file());
        let mut primes_text_file: File = File::open(primes_text_file_path).unwrap();
        let mut primes_text_content: String = String::new();
        let bytes_read = primes_text_file.read_to_string(&mut primes_text_content).unwrap();
        println!("Read {} bytes from prime table containing first million primes (primes1.txt)", bytes_read);
        return PrimeTableReader::read_primes(primes_text_content);
    }

    fn read_primes(primes_text_content: String) -> Vec<ZPlus> {
        let mut primes: Vec<ZPlus> = Vec::new();
        for line in primes_text_content.lines() {
            // Skip over empty lines
            if line.trim().is_empty() || line.contains("Primes") {
                continue;
            }

            for entry in line.split_whitespace() {
                let i_optional = entry.parse();
                if i_optional.is_ok() {
                    primes.push(i_optional.unwrap());
                } else {
                    panic!("Couldn't parse number from entry {}", entry)
                }
            }
        }

        return primes;
    }
}

#[cfg(test)]
mod prime_table_tests {
    use super::*;

    #[test]
    fn test_read_primes_to_one_million() {
        let reader: PrimeTableReader = PrimeTableReader::new(Path::new("prime_tables"));
        let primes: Vec<u64> = reader.first_million_primes();
        assert_eq!(1000000, primes.len())
    }
}