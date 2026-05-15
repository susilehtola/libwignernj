Name:           libwignernj
Version:        0.5.0
Release:        %autorelease
Summary:        Exact Wignernj 3j/6j/9j symbols and related coefficients via prime factorization
License:        BSD-3-Clause
URL:            https://github.com/susilehtola/libwignernj
Source0:        https://github.com/susilehtola/libwignernj/archive/v%{version}/%{name}-%{version}.tar.gz

BuildRequires:  cmake
BuildRequires:  gcc
BuildRequires:  gcc-c++
BuildRequires:  gcc-gfortran
BuildRequires:  libquadmath-devel
BuildRequires:  mpfr-devel
BuildRequires:  ninja-build
BuildRequires:  pkgconfig
BuildRequires:  python3-devel
# %%check
BuildRequires:  python3-pytest

%global _description %{expand:
libwignernj is a C99 library for the exact evaluation of Wignernj 3j, 6j and
9j symbols, Clebsch-Gordan coefficients, Racah W coefficients, Fano X
coefficients, and Gaunt coefficients (over both complex and real spherical
harmonics) using the prime-factorization scheme of Johansson and Forssen.
All intermediate arithmetic is exact integer arithmetic; floating-point
conversion happens only at the final step, so the returned values are
correct to the last representable bit of the chosen output type: float,
double, long double, IEEE 754 binary128 (via libquadmath), or, through
the optional GNU MPFR back-end, arbitrary precision.  Bindings are
provided for C, C++11 (header-only), Fortran 90 (iso_c_binding) and
Python 3.}

%description %_description

%package devel
Summary:        Development files for libwignernj
Requires:       libwignernj%{?_isa} = %{version}-%{release}
Requires:       mpfr-devel%{?_isa}
Requires:       libquadmath-devel%{?_isa}
Requires:       libstdc++-devel%{?_isa}

%description devel
C and C++ headers, pkg-config file, and CMake package configuration for
building applications that link against libwignernj.  Includes the
header-only C++11 wrapper (wignernj.hpp), the optional MPFR
arbitrary-precision header (wignernj_mpfr.h), and the optional libquadmath
IEEE 754 binary128 header (wignernj_quadmath.h).

%package fortran
Summary:        Fortran 90 bindings for libwignernj
Requires:       libwignernj%{?_isa} = %{version}-%{release}

%description fortran
Fortran 90 binding for libwignernj using the iso_c_binding intrinsic
module.  Provides real-valued wrappers (w3j, w6j, w9j, wcg, wracahw,
wfanox, wgaunt, wgaunt_real) as well as raw bind(c) interfaces in
single, double, long-double, and IEEE 754 binary128 (via libquadmath)
precisions.

%package fortran-devel
Summary:        Development files for the Fortran bindings to libwignernj
Requires:       libwignernj-fortran%{?_isa} = %{version}-%{release}
Requires:       libwignernj-devel%{?_isa} = %{version}-%{release}
Requires:       gcc-gfortran

%description fortran-devel
Fortran module (.mod) files needed to compile programs against the
Fortran 90 bindings of libwignernj.

%package -n python3-wignernj
Summary:        Python 3 bindings for libwignernj
%{?python_provide:%python_provide python3-wignernj}

%description -n python3-wignernj
Python 3 binding for libwignernj as a CPython extension dynamically
linked to libwignernj.so. Functions accept integer, float, or
fractions.Fraction arguments and dispatch to the appropriate C
precision through the precision='float'|'double'|'longdouble' keyword.

%package doc
Summary:        Documentation for libwignernj
BuildArch:      noarch

%description doc
The reference manual (Markdown) for libwignernj.

%prep
%setup -q

%build
%cmake \
    -DBUILD_FORTRAN:BOOL=ON \
    -DBUILD_MPFR:BOOL=ON \
    -DBUILD_QUADMATH:BOOL=ON \
    -DBUILD_PYTHON:BOOL=ON \
    -DBUILD_TESTS:BOOL=ON \
    -DBUILD_CXX_TESTS:BOOL=ON \
    -DBUILD_EXAMPLES:BOOL=ON \
    -GNinja
%cmake_build

%install
%cmake_install
%py_byte_compile %{python3} %{buildroot}%{python3_sitearch}/wignernj
mkdir -p %{buildroot}%{fmoddir}
mv %{buildroot}%{_includedir}/wignernj/fortran/*.mod %{buildroot}%{fmoddir}/

%check
%ctest
# Run the more comprehensive Python test suite, pointing pytest at the
# Python module and shared library that %install just placed in the
# buildroot (the system locations are not yet populated).
export PYTHONPATH=%{buildroot}%{python3_sitearch}
export LD_LIBRARY_PATH=%{buildroot}%{_libdir}
%pytest tests/python/

%files -n libwignernj
%license LICENSE
%doc README.md
%{_libdir}/libwignernj.so.0
%{_libdir}/libwignernj.so.0.*

%files devel
%{_includedir}/wignernj.h
%{_includedir}/wignernj.hpp
%{_includedir}/wignernj_mpfr.h
%{_includedir}/wignernj_quadmath.h
%{_libdir}/libwignernj.so
%{_libdir}/libwignernj_f03.so
%{_libdir}/cmake/wignernj/
%{_libdir}/pkgconfig/libwignernj.pc
%{_libdir}/pkgconfig/libwignernj_f03.pc

%files fortran
%license LICENSE
%{_libdir}/libwignernj_f03.so.0
%{_libdir}/libwignernj_f03.so.0.*


%files -n python3-wignernj
%license LICENSE
%doc README.md
%{python3_sitearch}/wignernj/

%files doc
%license LICENSE
%doc docs/reference.md
%doc docs/optimization_notes.md

%changelog
%autochangelog
