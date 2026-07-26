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

    if std::env::var_os("CARGO_FEATURE_LIBGAP").is_some() {
        let root = std::env::var("LIBGAP_ROOT")
            .expect("the libgap feature requires LIBGAP_ROOT to point to a built GAP checkout");
        let header = std::path::Path::new(&root).join("src/libgap-api.h");
        let library = std::path::Path::new(&root).join(if cfg!(target_os = "macos") {
            "libgap.dylib"
        } else {
            "libgap.so"
        });
        assert!(header.is_file(), "missing {}", header.display());
        assert!(library.is_file(), "missing {}", library.display());

        println!("cargo:rerun-if-env-changed=LIBGAP_ROOT");
        println!("cargo:rerun-if-changed=include/libgap_cyclotom.h");
        println!("cargo:rerun-if-changed=c/libgap_cyclotom.c");
        println!("cargo:rustc-link-search=native={root}");
        println!("cargo:rustc-link-lib=dylib=gap");
        println!("cargo:rustc-link-arg=-Wl,-rpath,{root}");

        cc::Build::new()
            .file("c/libgap_cyclotom.c")
            .include("include")
            .include(std::path::Path::new(&root).join("src"))
            .include(std::path::Path::new(&root).join("build"))
            .warnings(true)
            .extra_warnings(true)
            .flag_if_supported("-std=c11")
            .compile("libgap_cyclotom_bridge");
    }
}
