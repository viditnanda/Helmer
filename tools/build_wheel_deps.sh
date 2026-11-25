#!/bin/sh
# Build Helmer's C dependencies from source into $PREFIX (default /usr/local), for binary wheels:
# FLINT and msolve always, and GMP and MPFR when BUILD_GMP=1 (manylinux_2_28 ships GMP 6.1 and MPFR 3.1,
# too old for FLINT 3). Used by .github/workflows/wheels.yml.
#
# Portability. Wheels run on CPUs other than the one that built them, so nothing may use the build
# machine's instruction set:
#   * GMP is built "fat" (it picks its assembly at run time);
#   * FLINT 3.0 does not add -march=native;
#   * msolve's configure turns on every SIMD extension the build machine has (AVX2, AVX-512, BMI2, ...).
#     Its detection results are preset here so that it uses at most SSE4.2 (x86-64-v2, which numpy itself
#     now requires) and none of AVX, AVX-512, BMI, ADX, ABM, RDRND, SHA, AES.
set -eu
PREFIX=${PREFIX:-/usr/local}
JOBS=${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)}
GMP_VERSION=${GMP_VERSION:-6.3.0}
MPFR_VERSION=${MPFR_VERSION:-4.2.1}
FLINT_VERSION=${FLINT_VERSION:-3.0.1}     # the version Helmer is tested with (Ubuntu 24.04); msolve 0.6.5 needs a 3.0-era FLINT
MSOLVE_VERSION=${MSOLVE_VERSION:-0.6.5}
WORK=${WORK:-$(mktemp -d)}
cd "$WORK"
export CPPFLAGS="-I$PREFIX/include ${CPPFLAGS:-}" LDFLAGS="-L$PREFIX/lib ${LDFLAGS:-}"

if [ "${BUILD_GMP:-0}" = 1 ] && [ ! -f "$PREFIX/include/gmp.h" ]; then
  curl -fsSL -o gmp.tar.xz "https://ftp.gnu.org/gnu/gmp/gmp-$GMP_VERSION.tar.xz"
  tar xJf gmp.tar.xz && cd "gmp-$GMP_VERSION"
  ./configure --prefix="$PREFIX" --enable-fat --disable-static
  make -j"$JOBS" && make install && cd ..
  curl -fsSL -o mpfr.tar.xz "https://ftp.gnu.org/gnu/mpfr/mpfr-$MPFR_VERSION.tar.xz"
  tar xJf mpfr.tar.xz && cd "mpfr-$MPFR_VERSION"
  ./configure --prefix="$PREFIX" --with-gmp="$PREFIX" --disable-static
  make -j"$JOBS" && make install && cd ..
fi

if [ ! -f "$PREFIX/include/flint/flint.h" ]; then
  curl -fsSL -o flint.tar.gz "https://github.com/flintlib/flint/releases/download/v$FLINT_VERSION/flint-$FLINT_VERSION.tar.gz"
  tar xzf flint.tar.gz && cd "flint-$FLINT_VERSION"
  ./configure --prefix="$PREFIX" --disable-static ${GMP_PREFIX:+--with-gmp="$GMP_PREFIX" --with-mpfr="$GMP_PREFIX"}
  make -j"$JOBS" && make install && cd ..
fi

if [ ! -f "$PREFIX/include/msolve/neogb/libneogb.h" ]; then
  git clone --depth 1 --branch "v$MSOLVE_VERSION" https://github.com/algebraic-solving/msolve
  cd msolve
  ./autogen.sh
  ./configure --prefix="$PREFIX" --disable-static ${MSOLVE_CONFIGURE_EXTRA:-} \
      ax_cv_have_avx_os_support_ext=no ax_cv_have_avx512_os_support_ext=no \
      ax_cv_have_bmi1_cpu_ext=no ax_cv_have_bmi2_cpu_ext=no ax_cv_have_adx_cpu_ext=no \
      ax_cv_have_abm_cpu_ext=no ax_cv_have_rdrnd_cpu_ext=no ax_cv_have_prefetchwt1_cpu_ext=no \
      ax_cv_have_sha_cpu_ext=no ax_cv_have_aes_cpu_ext=no ax_cv_have_sse4a_cpu_ext=no
  if grep -E -- '-mavx|-mbmi|-madx|-march=native' Makefile src/*/Makefile; then
    echo "msolve would be built with CPU-specific flags" >&2; exit 1
  fi
  make -j"$JOBS" && make install && cd ..
fi
# Debug information makes the bundled libraries several times larger; wheels do not need it.
if command -v strip >/dev/null && [ "$(uname)" = Linux ]; then
  for f in "$PREFIX"/lib/lib*.so*; do [ -f "$f" ] && [ ! -L "$f" ] && strip --strip-debug "$f" || true; done
fi
echo "dependencies installed in $PREFIX"
