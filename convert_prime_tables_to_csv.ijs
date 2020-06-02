#!/usr/local/bin/ijconsole
load 'tables/csv'

firstMillionPrimesFile =: < 'prime_tables/primes1.txt'
readfile =: 1!:1
firstMillionPrimesText =: readfile firstMillionPrimesFile
newLineIndices =: I.@:(+./)@:((CR,LF)&(=/))
dropNewLines =: (<:@:[ { newLineIndices@:]) }. ]
firstMillionPrimesTextWithoutHeading =: 3 dropNewLines firstMillionPrimesText
firstMillionPrimesTextFlattened =: ' ' ([`(newLineIndices@:])`])} firstMillionPrimesTextWithoutHeading
firstMillionPrimesSpacedBoxes =: ' ' splitstring firstMillionPrimesTextFlattened
removeEmptyBoxes =: >@:(I.@:-.@:(*./"1)@:(' '&=)@:> { ])
firstMillionPrimesNumbers =: removeEmptyBoxes firstMillionPrimesSpacedBoxes
stringToNumber =: ".
firstMillionPrimes =: stringToNumber firstMillionPrimesNumbers
outputCsvFilename =: 'prime_tables/primes1.csv'
firstMillionPrimes writecsv outputCsvFilename
echo 'Wrote ', (": # firstMillionPrimes), ' primes to ',outputCsvFilename
exit 0