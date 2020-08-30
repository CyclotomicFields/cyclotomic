N := 100;;

G := CyclicGroup(N);;
ccs := ConjugacyClasses(G);;
sizes := List(ccs, Size);;
irr_chars := List(Irr(G), List);;
Print("sizes=", Concatenation(List(sizes, s -> Concatenation(" ", String(s)))), "\n");;
Print("num_chars=", Length(irr_chars), "\n");;
for char in irr_chars do
    Print(String(char), ";\n");;
od;;

# random character generation
rs := RandomSource(IsMersenneTwister, NanosecondsSinceEpoch());;
char := Sum([1..Length(ccs)], i -> Random(rs, [-20..20]) * irr_chars[i]);;
Print("random_char=", String(char), ";\n");;

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

PrintTo("*errout*", "GAP calculated prod: ", prods, "\n");
ms_elapsed := Float((end_time-start_time)/1000000);;
PrintTo("*errout*", "GAP time elapsed: ", ms_elapsed, "ms\n");

AppendTo("gap_results", N, ",", ms_elapsed, "\n");;
