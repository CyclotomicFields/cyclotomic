use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

type ZPlus = u64;

pub struct PrimeTableReader<'a> {
    directory_path: &'a Path,
}

impl<'a> PrimeTableReader<'a> {
    pub fn new(directory_path: &Path) -> PrimeTableReader {
        assert!(directory_path.is_dir());
        let first_million_primes_txt_path = directory_path.join("primes1.txt");
        let first_million_primes_csv_path = directory_path.join("primes1.csv");
        assert!(first_million_primes_txt_path.is_file() || first_million_primes_csv_path.is_file());
        return PrimeTableReader { directory_path };
    }

    pub fn first_million_from_file() -> Option<PrimeTableReader<'a>> {
        let prime_tables_path = Path::new("prime_tables");
        let first_million_primes_txt_path = prime_tables_path.join("primes1.txt");
        let first_million_primes_csv_path = prime_tables_path.join("primes1.csv");
        return if prime_tables_path.is_dir()
            && (first_million_primes_txt_path.is_file() || first_million_primes_csv_path.is_file())
        {
            Some(PrimeTableReader::new(prime_tables_path))
        } else {
            None
        };
    }

    pub fn first_million_primes(&self) -> Vec<ZPlus> {
        let primes_csv_file_path: PathBuf = self.directory_path.join("primes1.csv");
        let primes_txt_file_path: PathBuf = self.directory_path.join("primes1.txt");
        return if primes_csv_file_path.is_file() {
            let mut primes_text_file: File = File::open(primes_csv_file_path).unwrap();
            let mut primes_text_content: String = String::new();
            let bytes_read = primes_text_file
                .read_to_string(&mut primes_text_content)
                .unwrap();
            println!(
                "Read {} bytes from prime table containing first million primes (primes1.csv)",
                bytes_read
            );
            PrimeTableReader::read_primes_from_csv(primes_text_content)
        } else if primes_txt_file_path.is_file() {
            let mut primes_text_file: File = File::open(primes_txt_file_path).unwrap();
            let mut primes_text_content: String = String::new();
            let bytes_read = primes_text_file
                .read_to_string(&mut primes_text_content)
                .unwrap();
            println!(
                "Read {} bytes from prime table containing first million primes (primes1.txt)",
                bytes_read
            );
            PrimeTableReader::read_primes_from_table(primes_text_content)
        } else {
            panic!(
                "Couldn't find either {:?} or {:?} inside directory {}",
                primes_csv_file_path,
                primes_txt_file_path,
                self.directory_path.display()
            )
        };
    }

    fn read_primes_from_table(primes_txt_content: String) -> Vec<ZPlus> {
        let mut primes: Vec<ZPlus> = Vec::new();
        for line in primes_txt_content.lines() {
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

    fn read_primes_from_csv(primes_csv_content: String) -> Vec<ZPlus> {
        let mut primes: Vec<ZPlus> = Vec::new();
        for entry in primes_csv_content.split(",") {
            let i_optional = entry.trim().parse();
            if i_optional.is_ok() {
                primes.push(i_optional.unwrap());
            } else {
                panic!("Couldn't parse number from entry {}", entry)
            }
        }
        return primes;
    }
}

// Runs only if primes1.txt has been downloaded.
#[cfg(test)]
mod prime_table_tests {
    use super::*;

    #[test]
    fn test_read_primes_to_one_million() {
        if let Some(prime_table_reader) = PrimeTableReader::first_million_from_file() {
            assert_eq!(1000000, prime_table_reader.first_million_primes().len())
        }
    }
}
