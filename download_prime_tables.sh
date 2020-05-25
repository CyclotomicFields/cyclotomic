first_fifty_million_primes="https://primes.utm.edu/lists/small/millions/primes1.zip"
curl $first_fifty_million_primes --output download.zip
unzip -o download.zip
mv primes1.txt prime_tables
rm download.zip
