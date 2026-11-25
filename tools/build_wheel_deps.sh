#!/bin/sh
# Build FLINT and msolve from source into $PREFIX (default /usr/local), for binary wheels.
# GMP and MPFR come from the system (manylinux: dnf; macOS: Homebrew). Used by .github/workflows/wheels.yml.
#
# msolve's configure enables every SIMD extension of the build machine (AVX2, AVX-512, BMI2, ...). A wheel
# built that way can crash with an illegal instruction on an older CPU, so CPU-specific flags are fixed here:
# MSOLVE_SIMD (default empty: the baseline of the architecture, SSE2 on x86-64).
set -eu
PREFIX=${PREFIX:-/usr/local}
JOBS=${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)}
FLINT_VERSION=${FLINT_VERSION:-3.1.2}
MSOLVE_VERSION=${MSOLVE_VERSION:-0.6.5}
MSOLVE_SIMD=${MSOLVE_SIMD:-}
WORK=${WORK:-$(mktemp -d)}
cd "$WORK"

if [ ! -f "$PREFIX/include/flint/flint.h" ]; then
  curl -fsSL -o flint.tar.gz "https://github.com/flintlib/flint/releases/download/v$FLINT_VERSION/flint-$FLINT_VERSION.tar.gz"
  tar xzf flint.tar.gz
  cd "flint-$FLINT_VERSION"
  # no -march=native: a portable library
  ./configure --prefix="$PREFIX" --disable-static ${FLINT_CONFIGURE_EXTRA:-}
  make -j"$JOBS"
  make install
  cd ..
fi

if [ ! -f "$PREFIX/include/msolve/neogb/libneogb.h" ]; then
  git clone --depth 1 --branch "v$MSOLVE_VERSION" https://github.com/algebraic-solving/msolve
  cd msolve
  ./autogen.sh
  ./configure --prefix="$PREFIX" --disable-static \
      CPPFLAGS="-I$PREFIX/include ${CPPFLAGS:-}" LDFLAGS="-L$PREFIX/lib ${LDFLAGS:-}" ${MSOLVE_CONFIGURE_EXTRA:-}
  make -j"$JOBS" SIMD_FLAGS="$MSOLVE_SIMD" CPUEXT_FLAGS=""
  make install
  cd ..
fi
echo "FLINT $FLINT_VERSION and msolve $MSOLVE_VERSION installed in $PREFIX"
