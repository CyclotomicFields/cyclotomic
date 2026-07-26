fn main() {
    println!("cargo:rerun-if-changed=include/gap_cyclotom.h");
    println!("cargo:rerun-if-changed=c/gap_cyclotom.c");

    cc::Build::new()
        .file("c/gap_cyclotom.c")
        .include("include")
        .warnings(true)
        .extra_warnings(true)
        .flag_if_supported("-std=c11")
        .compile("gap_cyclotom_reference");
}
