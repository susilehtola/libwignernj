# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
#
# Build script for the comparative benchmark in this directory.  Three
# external dependencies are required:
#
#   - libwignernj          (this repository; built and either installed
#                           or pointed to via WIGNERNJ_DIR)
#   - WIGXJPF >= 1.13      (https://fy.chalmers.se/subatom/wigxjpf/);
#                           pass WIGXJPF_DIR=/path/to/wigxjpf
#   - GNU Scientific Library (system package, located via pkg-config)
#
# Usage:
#
#   make WIGXJPF_DIR=/path/to/wigxjpf-1.13 \
#        WIGNERNJ_DIR=/path/to/wignernj/build
#
# Or, with both libraries installed system-wide so that pkg-config and
# the default include/lib search paths resolve them:
#
#   make
#
# The comparative benchmark in the published paper used the default
# Fedora optimization flags, which on a Fedora system can be queried
# with `rpm -E %optflags`.  The same flags should be used to build all
# three libraries (libwignernj, WIGXJPF, and GSL) for a fair
# comparison.  Pass CFLAGS=... on the command line to override.
#
# When this Makefile runs on a Fedora-derived system with rpm
# installed, the default CFLAGS are taken directly from
# `rpm -E %optflags`; otherwise an explicit literal copy is used as
# a fallback.

CFLAGS ?= $(shell rpm -E %optflags 2>/dev/null \
          || echo "-O2 -flto=auto -ffat-lto-objects -fexceptions -g -grecord-gcc-switches \
                   -pipe -Wall -Wno-complain-wrong-lang -Werror=format-security \
                   -Wp,-U_FORTIFY_SOURCE,-D_FORTIFY_SOURCE=3 \
                   -Wp,-D_GLIBCXX_ASSERTIONS \
                   -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1 \
                   -fstack-protector-strong \
                   -specs=/usr/lib/rpm/redhat/redhat-annobin-cc1 \
                   -m64 -march=x86-64 -mtune=generic \
                   -fasynchronous-unwind-tables -fstack-clash-protection -fcf-protection \
                   -mtls-dialect=gnu2 -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer")

# libwignernj: prefer pkg-config, fall back to WIGNERNJ_DIR.
ifdef WIGNERNJ_DIR
  WIGNERNJ_CFLAGS := -I$(WIGNERNJ_DIR)/../include -I../include
  WIGNERNJ_LDFLAGS := -L$(WIGNERNJ_DIR) -lwignernj
else
  WIGNERNJ_CFLAGS := $(shell pkg-config --cflags libwignernj)
  WIGNERNJ_LDFLAGS := $(shell pkg-config --libs libwignernj)
endif

# WIGXJPF: WIGXJPF_DIR points to the unpacked source where ./inc/ and
# ./lib/ live (i.e. wigxjpf-1.13/).
ifndef WIGXJPF_DIR
  $(error Please pass WIGXJPF_DIR=/path/to/wigxjpf-1.13)
endif
WIGXJPF_CFLAGS  := -I$(WIGXJPF_DIR)/inc
WIGXJPF_LDFLAGS := $(WIGXJPF_DIR)/lib/libwigxjpf.a

# GSL via pkg-config.
GSL_CFLAGS  := $(shell pkg-config --cflags gsl)
GSL_LDFLAGS := $(shell pkg-config --libs   gsl)

CC ?= cc

bench_compare: bench_compare.c
	$(CC) $(CFLAGS) $(WIGNERNJ_CFLAGS) $(WIGXJPF_CFLAGS) $(GSL_CFLAGS) \
	      bench_compare.c \
	      $(WIGNERNJ_LDFLAGS) $(WIGXJPF_LDFLAGS) $(GSL_LDFLAGS) -lm \
	      -o $@

run: bench_compare
	./bench_compare

clean:
	rm -f bench_compare

.PHONY: run clean
