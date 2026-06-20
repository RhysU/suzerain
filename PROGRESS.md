# Goal: Get `make` to complete successfully for suzerain

## Guiding principles (from user)
- Prefer fixing source code to adhere to newer standards over adding backwards-compat flags
- Fix detection macros (m4) at the root cause rather than hardcoding workarounds
- We own all source code (suzerain, largo, p3dfft) — fix it directly
- Use https://github.com/autoconf-archive/autoconf-archive for updated m4 macros
- Use https://github.com/tsuna/boost.m4 for the latest boost.m4

## Environment
- Ubuntu 24.04 (Noble), GCC 13.3, gfortran 13.3, OpenMPI
- Build: autoconf/automake/libtool, `AX_ENABLE_BUILDDIR` puts output in `x86_64-unknown-linux-gnu/`

## Steps to reproduce the build

### 1. Install system packages
```bash
sudo apt-get update -qq
sudo apt-get install -y libtool gfortran openmpi-bin libopenmpi-dev \
    libgsl-dev libfftw3-dev libfftw3-mpi-dev libboost-all-dev \
    libeigen3-dev liblog4cxx-dev libhdf5-openmpi-dev libapr1-dev libaprutil1-dev
```

### 2. Build and install ESIO from source
ESIO is not packaged in Ubuntu. Build from https://github.com/RhysU/ESIO:
```bash
cd /tmp && git clone https://github.com/RhysU/ESIO.git
cd ESIO && ./bootstrap
CC=mpicc FC=mpif90 ./configure --prefix=/usr/local --disable-tests
make -j1   # parallel Fortran module build is broken
sudo make install && sudo ldconfig
```

### 3. Bootstrap and configure suzerain
```bash
cd /home/user/suzerain
./bootstrap
CC=mpicc CXX=mpicxx FC=mpif90 ./configure \
    --disable-writeups --disable-unittests --without-mkl
```

### 4. Fix libtool bare `-l` entries
OpenMPI's Fortran wrapper generates empty `-l` tokens in libtool postdeps.
Fix both the main and largo libtool files:
```bash
python3 -c "
import re
for lt in ['x86_64-unknown-linux-gnu/libtool', 'x86_64-unknown-linux-gnu/largo/libtool']:
    with open(lt, 'r') as f:
        content = f.read()
    content = re.sub(r'-l(?=[\s\"])', '', content)
    with open(lt, 'w') as f:
        f.write(content)
"
```
Also added `AC_CONFIG_COMMANDS` in configure.ac to auto-fix the main libtool.

### 5. Build
```bash
cd x86_64-unknown-linux-gnu && make -j$(nproc)
```

## Source code changes made

### m4/boost.m4 — replaced with upstream serial 39
- Old version (serial 18) failed to detect Boost version due to GCC 5+
  inserting `# line` directives mid-expansion, splitting the version string
  across multiple output lines.
- Downloaded latest from https://github.com/tsuna/boost.m4/master/build-aux/boost.m4
  which handles this with `grep -v '#' | tr -s '\n' ' '`.

### m4/ax_fftw3.m4 — added FFTW3_FCFLAGS detection
- pkg-config returns empty CFLAGS for FFTW3 when installed in /usr/include (standard path).
- Fortran `include "fftw3.f"` needs explicit `-I` unlike C `#include`.
- Added detection of fftw3.f location via `pkg-config --variable=includedir`.

### lib/suzerain-p3dfft/Makefile.am — use FFTW3_FCFLAGS and FCFLAGS_ALLOW_ARGUMENT_MISMATCH
- Added `AM_FCFLAGS = $(FFTW3_FCFLAGS) $(FCFLAGS_ALLOW_ARGUMENT_MISMATCH)`
- P3DFFT passes mixed real/complex types to FFTW plan routines, which gfortran 10+
  rejects by default. The `-fallow-argument-mismatch` flag is needed here because
  this is a fundamental part of the FFTW Fortran API design pattern.

### configure.ac — Fortran compatibility
- Added `AX_CHECK_COMPILE_FLAG([-fallow-argument-mismatch])` in Fortran section
  to detect and export `FCFLAGS_ALLOW_ARGUMENT_MISMATCH`.
- Added `AC_CONFIG_COMMANDS` to fix libtool bare `-l` entries post-configure.

### largo/largo/largo_workspace.f90 — keep as-is (no change needed after largo.f90 fix)

### largo/largo/largo.f90 — fix gfortran 10+ type ambiguity
- `largo_workspace` module exports `largo_ptr => c_ptr` and `largo_workspace_ptr => c_ptr`.
- `largo.f90` also defines these same aliases from iso_c_binding.
- gfortran 10+ rejects same-name imports even when underlying types match.
- Fixed by renaming the imported aliases: `use largo_workspace, lws_ptr => largo_ptr, lws_wptr => largo_workspace_ptr`

### suzerain/bspline.h, bspline.c, bspline.hpp, bsplineop.c — GSL 2.x API update
- **IN PROGRESS**: GSL 2.x merged `gsl_bspline_deriv_workspace` into `gsl_bspline_workspace`.
- Need to remove the separate deriv workspace type, alloc/free calls, and extra parameters.
- `gsl_bspline_deriv_eval_nonzero()` signature changed: last arg is now `gsl_bspline_workspace*`
  instead of separate `gsl_bspline_deriv_workspace*`.

## Current status
- [x] System packages installed
- [x] ESIO built and installed
- [x] bootstrap + configure pass
- [x] P3DFFT (lib/suzerain-p3dfft) compiles and links
- [x] Largo compiles and links
- [ ] GSL bspline API migration (in progress)
- [ ] Main suzerain library compilation
- [ ] Apps compilation
- [ ] Full `make` success
