N := 54;;

G := CyclicGroup(N);;
ccs := ConjugacyClasses(G);;
sizes := List(ccs, Size);;
irr_chars := List(Irr(G), List);;
PrintFormattedString(Concatenation("sizes=", Concatenation(List(sizes, s -> Concatenation(" ", String(s)))), "\n"));;
PrintFormattedString(Concatenation("num_chars=", String(Length(irr_chars)), "\n"));;
for char in irr_chars do
    PrintFormattedString(Concatenation(String(char), ";\n"));;
od;;

# random character generation
rs := RandomSource(IsMersenneTwister, NanosecondsSinceEpoch());;
char := Sum([1..Length(ccs)], i -> Random(rs, [1..2]) * irr_chars[i]);;

# just pick a single one for max sparsity
#char := Random(rs, irr_chars);;

PrintFormattedString(Concatenation("random_char=", String(char), ";\n"));;

# character inner product
prod := function(char1, char2, sizes)
    local i, n;;
    n := Length(char1);;
    return Sum([1..n], i -> sizes[i] * char1[i] * ComplexConjugate(char2[i]));;
end;;

# run benchmark
start_time := NanosecondsSinceEpoch();;

prods := List(irr_chars, irr_char -> prod(irr_char, char, sizes));;

end_time := NanosecondsSinceEpoch();;

#PrintTo("*errout*", "GAP calculated prod: ", prods, "\n");
ns_elapsed := (end_time-start_time);;
PrintTo("*errout*", "GAP prods:: ", prods, "\n");
PrintTo("*errout*", "GAP time elapsed: ", ns_elapsed, " ns\n");

AppendTo("gap_results", N, ",", ns_elapsed/1000000, "\n");;
