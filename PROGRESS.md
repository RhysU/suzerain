# Goal: Get `make` to complete successfully for suzerain

## Guiding principles (from user)
- Prefer fixing source code to adhere to newer standards over adding backwards-compat flags
- Fix detection macros (m4) at the root cause rather than hardcoding workarounds
- We own all source code (suzerain, largo, p3dfft) — fix it directly
- Use https://github.com/autoconf-archive/autoconf-archive for updated m4 macros
- Use https://github.com/tsuna/boost.m4 for the latest boost.m4
- Prefer to install missing prerequisites rather than working around their absence
- Be wary of changes requiring cascading changes
- DO NOT include Claude or Anthropic in commit authorship

## Environment
- Ubuntu 24.04 (Noble), GCC 13.3, gfortran 13.3, OpenMPI
- Build: autoconf/automake/libtool, `AX_ENABLE_BUILDDIR` puts output in `x86_64-unknown-linux-gnu/`

## Steps to reproduce the build

### 1. Install system packages
```bash
sudo apt-get update -qq
sudo apt-get install -y libtool gfortran openmpi-bin libopenmpi-dev \
    libgsl-dev libfftw3-dev libfftw3-mpi-dev libboost-all-dev \
    libeigen3-dev liblog4cxx-dev libhdf5-openmpi-dev libapr1-dev libaprutil1-dev \
    libmkl-dev
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
    --disable-writeups --disable-unittests \
    --with-mkl-include=/usr/include/mkl --with-mkl-libdir=/usr/lib/x86_64-linux-gnu
```

### 4. Fix libtool bare `-l` entries and dep linking
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
Also fix `link_all_deplibs=no` to `link_all_deplibs=unknown` in libtool
to ensure bare `-l` dependency libs from `.la` files are passed through:
```bash
sed -i 's/^link_all_deplibs=no$/link_all_deplibs=unknown/' x86_64-unknown-linux-gnu/libtool
```

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

### m4/acx_mkl.m4 — use mkl_get_version instead of MKLGetVersion
- Ubuntu's MKL 2020.4 (libmkl-dev) does not export `MKLGetVersion`.
- Replaced with `mkl_get_version` which is the current MKL API function.
- Important to check for MKL specifically (not just dgemm which could come from other BLAS).

### lib/suzerain-p3dfft/Makefile.am — use FFTW3_FCFLAGS and FCFLAGS_ALLOW_ARGUMENT_MISMATCH
- Added `AM_FCFLAGS = $(FFTW3_FCFLAGS) $(FCFLAGS_ALLOW_ARGUMENT_MISMATCH)`
- P3DFFT passes mixed real/complex types to FFTW plan routines, which gfortran 10+
  rejects by default. The `-fallow-argument-mismatch` flag is needed here because
  this is a fundamental part of the FFTW Fortran API design pattern.

### configure.ac — Fortran compatibility
- Added `AX_CHECK_COMPILE_FLAG([-fallow-argument-mismatch])` in Fortran section
  to detect and export `FCFLAGS_ALLOW_ARGUMENT_MISMATCH`.
- Added `AC_CONFIG_COMMANDS` to fix libtool bare `-l` entries post-configure.

### largo/largo/largo.f90 — fix gfortran 10+ type ambiguity
- `largo_workspace` module exports `largo_ptr => c_ptr` and `largo_workspace_ptr => c_ptr`.
- `largo.f90` also defines these same aliases from iso_c_binding.
- gfortran 10+ rejects same-name imports even when underlying types match.
- Fixed by renaming the imported aliases.

### suzerain/bspline.h, bspline.c, bspline.hpp — GSL 2.x bspline API migration
- GSL 2.x merged `gsl_bspline_deriv_workspace` into `gsl_bspline_workspace`.
- Removed the separate `dw` parameter from all bspline functions.
- `gsl_bspline_deriv_eval_nonzero()` last arg is now `gsl_bspline_workspace*`.

### suzerain/bsplineop.h, bsplineop.c — GSL 2.x bspline API migration
- Removed `gsl_bspline_deriv_workspace *dw` from `suzerain_bsplineop_workspace` struct.
- Removed alloc/free of the deriv workspace.
- Updated all calls to pass `w` instead of `w, dw`.

### suzerain/bl.h, bl.c — GSL 2.x bspline API migration
- Removed `gsl_bspline_deriv_workspace *dw` parameter from 4 public functions:
  `suzerain_bl_find_edge`, `suzerain_bl_find_edge99`,
  `suzerain_bl_compute_thicknesses`, `suzerain_bl_compute_thicknesses_baseflow`.

### suzerain/mpi.hpp, mpi.cpp — C++17 compatibility
- Removed `throw(std::logic_error)` dynamic exception specification from
  `ensure_mpi_initialized()` (illegal in C++17).

### suzerain/validation.hpp — C++17 compatibility
- Removed `throw(std::invalid_argument)` dynamic exception specifications
  from `ensure_positive`, `ensure_nonnegative`, and `ensure_bounded`.

### suzerain/complex.hpp — C++17 ADL ambiguity fix
- `suzerain::complex::real(const FPT&)` conflicts with `std::real(T)` in C++17.
- Removed `using std::real;` and `using std::imag;` directives.
- Added explicit `std::complex<FPT>` overloads returning lvalue references.
- Qualified all internal calls with `(real)(x)` / `(imag)(x)` idiom to suppress ADL.
- Future work: migrate to standard C++ real/imag (see GitHub issue #1).

### suzerain/utility.hpp — remove ambiguous EigenDenseType overloads
- `to_yxz` and `to_xzy` had both `RandomAccessContainer` and `EigenDenseType`
  template overloads that became ambiguous because modern Eigen defines
  `value_type = Scalar` in `DenseBase`.
- Removed the `EigenDenseType` overloads since the `RandomAccessContainer`
  versions work for both (Eigen types have `value_type`).

### suzerain/Makefile.am — isolate MKL include path
- Changed `libbspline_la_CPPFLAGS` to inherit from `libsuzerain_la_CPPFLAGS`
  instead of `libblasetal_la_CPPFLAGS` to prevent `-I/usr/include/mkl` from
  reaching code that includes GSL headers.
- Avoids MKL CBLAS / GSL CBLAS header conflict.

### suzerain/support/application_base.cpp — MKL cleanup and fixes
- Added `extern "C" void mkl_free_buffers(void);` declaration instead of
  including mkl.h (which conflicts with GSL's CBLAS headers).
- Simplified MKL buffer cleanup code (removed old `INTEL_MKL_VERSION` branch).
- Fixed ambiguous `to_yxz` call with `static_cast<int>(linear_nfields)`.

### suzerain/support/logging.cpp — log4cxx 1.1 API migration
- Removed `ObjectImpl` inheritance (removed in log4cxx 1.1).
- Updated `SubversiveASHEL` to use raw pointer signatures per log4cxx 1.1's
  `HierarchyEventListener` interface.
- Changed `new SubversiveASHEL()` to `std::make_shared<SubversiveASHEL>()`.
- MPI rank awareness preserved: rank > 0 intercepts appenders and redirects
  root logger appenders to rank-specific logger subtrees.

### apps/perfect/perfect.cpp — updated call sites
- Removed `b.dbw` from `suzerain_bl_compute_thicknesses` and
  `suzerain_bl_compute_thicknesses_baseflow` calls.

### apps/perfect/main_initial.cpp — disambiguate make_shared
- Qualified `boost::make_shared<manufactured_solution>(...)` to resolve
  ambiguity with `std::make_shared` found via ADL on `std::string` argument.

## Known issues / future work
- `boost::make_shared` vs `std::make_shared` ambiguity: the codebase uses
  `boost::shared_ptr` throughout (via `using boost::shared_ptr` in common.hpp).
  Cannot switch to `std::make_shared` without also migrating all shared_ptr usage.
- complex.hpp real/imag migration to standard C++ (GitHub issue #1)
- libtool `link_all_deplibs` fix is needed after each `configure` run

## Current status
- [x] System packages installed (including libmkl-dev)
- [x] ESIO built and installed
- [x] bootstrap + configure pass (with MKL)
- [x] P3DFFT (lib/suzerain-p3dfft) compiles and links
- [x] Largo compiles and links
- [x] GSL bspline API migration complete
- [x] Main suzerain library compilation
- [x] Apps compilation
- [x] **Full `make` success**
