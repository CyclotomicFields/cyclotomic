#!/usr/bin/env bash
set -euo pipefail

commit=d2134de71521c62512b8351c42ec16bfbac21744
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
default_root="$script_dir/../target/libgap-$commit"
gap_root=${1:-$default_root}

if [[ ! -d "$gap_root/.git" ]]; then
    git clone --filter=blob:none https://github.com/gap-system/gap.git "$gap_root"
fi

git -C "$gap_root" fetch origin "$commit"
git -C "$gap_root" checkout --detach "$commit"

configure_args=()
if command -v brew >/dev/null 2>&1; then
    gmp_prefix=$(brew --prefix gmp 2>/dev/null || true)
    if [[ -n "$gmp_prefix" ]]; then
        configure_args+=("--with-gmp=$gmp_prefix")
    fi
fi

(
    cd "$gap_root"
    ./autogen.sh
    ./configure "${configure_args[@]}"
    make -j"${CYCLOTOMIC_BUILD_JOBS:-4}"
    make bootstrap-pkg-minimal
)

if [[ $(uname -s) == Darwin ]]; then
    libgap_target=$(readlink "$gap_root/libgap.dylib")
    install_name_tool -id "@rpath/$libgap_target" "$gap_root/$libgap_target"
fi

printf 'export LIBGAP_ROOT=%q\n' "$gap_root"
