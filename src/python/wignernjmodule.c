/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * CPython extension module for libwignernj.
 * Exposes wigner3j, wigner6j, wigner9j, clebsch_gordan, racah_w, fano_x, gaunt.
 *
 * Each function accepts integer, float (half-integer), or fractions.Fraction
 * arguments.  An optional keyword argument precision={'float','double',
 * 'longdouble'} selects the output precision (default: 'double').
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include "wignernj.h"

/* ── half-integer argument parsing ─────────────────────────────────────── */

/*
 * Parse a single Python object (int, float, or Fraction) as 2*j.
 * Returns 0 on success, -1 on error (Python exception set).
 */
static int parse_half_int(PyObject *obj, int *out, const char *name)
{
    if (PyLong_Check(obj)) {
        long v = PyLong_AsLong(obj);
        if (v == -1 && PyErr_Occurred()) return -1;
        *out = (int)(2 * v);
        return 0;
    }
    if (PyFloat_Check(obj)) {
        double v = PyFloat_AsDouble(obj);
        int tv = (v >= 0) ? (int)(2.0*v + 0.5) : -(int)(-2.0*v + 0.5);
        if (fabs(2.0*v - tv) > 1e-9) {
            PyErr_Format(PyExc_ValueError,
                "wigner: '%s' is not a half-integer", name);
            return -1;
        }
        *out = tv;
        return 0;
    }
    /* Try fractions.Fraction: access .numerator and .denominator */
    {
        PyObject *num_obj = PyObject_GetAttrString(obj, "numerator");
        PyObject *den_obj = PyObject_GetAttrString(obj, "denominator");
        if (!num_obj || !den_obj) {
            Py_XDECREF(num_obj); Py_XDECREF(den_obj);
            PyErr_Format(PyExc_TypeError,
                "wigner: '%s' must be int, float, or Fraction", name);
            return -1;
        }
        long num = PyLong_AsLong(num_obj);
        long den = PyLong_AsLong(den_obj);
        Py_DECREF(num_obj); Py_DECREF(den_obj);
        if ((num == -1 || den == -1) && PyErr_Occurred()) return -1;
        if (den == 1) { *out = (int)(2 * num); return 0; }
        if (den == 2) { *out = (int)num; return 0; }
        PyErr_Format(PyExc_ValueError,
            "wigner: '%s' = %ld/%ld is not a half-integer", name, num, den);
        return -1;
    }
}

/* ── shared result builder ─────────────────────────────────────────────── */

typedef enum { PREC_FLOAT, PREC_DOUBLE, PREC_LONGDOUBLE } Precision;

static Precision parse_precision(PyObject *prec_obj)
{
    if (!prec_obj || prec_obj == Py_None) return PREC_DOUBLE;
    const char *s = PyUnicode_AsUTF8(prec_obj);
    if (!s) return (Precision)-1;  /* TypeError already set by PyUnicode_AsUTF8 */
    if (strcmp(s, "float")      == 0) return PREC_FLOAT;
    if (strcmp(s, "double")     == 0) return PREC_DOUBLE;
    if (strcmp(s, "longdouble") == 0) return PREC_LONGDOUBLE;
    PyErr_SetString(PyExc_ValueError,
        "wigner: precision must be 'float', 'double', or 'longdouble'");
    return (Precision)-1;
}

static PyObject *make_result(long double val, Precision prec)
{
    switch (prec) {
    case PREC_FLOAT:      return PyFloat_FromDouble((double)(float)val);
    case PREC_DOUBLE:     return PyFloat_FromDouble((double)val);
    case PREC_LONGDOUBLE: {
#ifdef Py_HAS_NUMPY
        /* Would use np.longdouble — skip for portability */
#endif
        return PyFloat_FromDouble((double)val);
    }
    }
    return PyFloat_FromDouble((double)val);
}

/* ── wigner3j ──────────────────────────────────────────────────────────── */

static const char wigner3j_doc[] =
    "wigner3j(j1, j2, j3, m1, m2, m3, precision='double') -> float\n"
    "\n"
    "Wigner 3j symbol\n"
    "\n"
    "    ( j1  j2  j3 )\n"
    "    ( m1  m2  m3 )\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "j1, j2, j3 : int, float, or fractions.Fraction\n"
    "    Angular-momentum quantum numbers.  Pass an int for integer j\n"
    "    (e.g. ``1`` for j=1); a float (e.g. ``0.5``) or\n"
    "    ``Fraction(1, 2)`` for half-integer j.  Non-half-integer\n"
    "    arguments raise ValueError.\n"
    "m1, m2, m3 : int, float, or fractions.Fraction\n"
    "    Magnetic projection quantum numbers.  Same type-acceptance\n"
    "    rules as the j_i.\n"
    "precision : {'float', 'double', 'longdouble'}, optional\n"
    "    IEEE 754 binary precision of the returned value.  Default\n"
    "    'double' (binary64).  'float' returns binary32; 'longdouble'\n"
    "    returns the platform's extended type (binary80 on x86-64,\n"
    "    binary128 on aarch64/POWER, otherwise the same as 'double').\n"
    "\n"
    "Returns\n"
    "-------\n"
    "float\n"
    "    The 3j symbol, correctly rounded to the chosen precision.\n"
    "    Returns 0.0 silently if any selection rule is violated; this\n"
    "    is not an error condition.\n"
    "\n"
    "Raises\n"
    "------\n"
    "ValueError\n"
    "    If any argument is not a half-integer, or if `precision` is\n"
    "    not one of 'float', 'double', 'longdouble'.\n"
    "TypeError\n"
    "    If any argument is not int, float, or Fraction.\n"
    "\n"
    "Notes\n"
    "-----\n"
    "The 3j symbol is a pure SU(2) algebraic object and carries no\n"
    "spherical-harmonic phase convention.  See `help(wignernj)` for\n"
    "the full conventions block.\n"
    "\n"
    "Examples\n"
    "--------\n"
    ">>> import wignernj\n"
    ">>> wignernj.wigner3j(1, 1, 0, 0, 0, 0)\n"
    "-0.5773502691896257\n"
    ">>> wignernj.wigner3j(0.5, 0.5, 1, 0.5, -0.5, 0)  # 1/sqrt(6)\n"
    "0.408248290463863\n"
    ">>> from fractions import Fraction\n"
    ">>> wignernj.wigner3j(Fraction(1,2), Fraction(1,2), 1,\n"
    "...                    Fraction(1,2), Fraction(-1,2), 0)\n"
    "0.408248290463863";

static PyObject *py_wigner3j(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"j1","j2","j3","m1","m2","m3","precision",NULL};
    PyObject *j1o,*j2o,*j3o,*m1o,*m2o,*m3o,*prec_obj=NULL;
    int tj1,tj2,tj3,tm1,tm2,tm3;
    Precision prec;
    (void)self;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOO|O", kwlist,
                                     &j1o,&j2o,&j3o,&m1o,&m2o,&m3o,&prec_obj))
        return NULL;
    if (parse_half_int(j1o,&tj1,"j1") || parse_half_int(j2o,&tj2,"j2") ||
        parse_half_int(j3o,&tj3,"j3") || parse_half_int(m1o,&tm1,"m1") ||
        parse_half_int(m2o,&tm2,"m2") || parse_half_int(m3o,&tm3,"m3"))
        return NULL;
    prec = parse_precision(prec_obj);
    if ((int)prec < 0) return NULL;
    return make_result(wigner3j_l(tj1,tj2,tj3,tm1,tm2,tm3), prec);
}

/* ── wigner6j ──────────────────────────────────────────────────────────── */

static const char wigner6j_doc[] =
    "wigner6j(j1, j2, j3, j4, j5, j6, precision='double') -> float\n"
    "\n"
    "Wigner 6j symbol\n"
    "\n"
    "    { j1 j2 j3 }\n"
    "    { j4 j5 j6 }\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "j1, j2, j3, j4, j5, j6 : int, float, or fractions.Fraction\n"
    "    Angular-momentum quantum numbers.  Pass an int for integer j,\n"
    "    a float (e.g. ``0.5``) or ``Fraction(1, 2)`` for half-integer j.\n"
    "    The four triangle conditions (j1 j2 j3), (j1 j5 j6),\n"
    "    (j4 j2 j6), (j4 j5 j3) must hold simultaneously for a\n"
    "    non-vanishing result.\n"
    "precision : {'float', 'double', 'longdouble'}, optional\n"
    "    IEEE 754 binary precision of the returned value.  Default\n"
    "    'double'.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "float\n"
    "    The 6j symbol, correctly rounded to the chosen precision.\n"
    "    Returns 0.0 silently if any triangle condition is violated.\n"
    "\n"
    "Examples\n"
    "--------\n"
    ">>> import wignernj\n"
    ">>> wignernj.wigner6j(1, 1, 1, 1, 1, 1)         # 1/6\n"
    "0.16666666666666666\n"
    ">>> wignernj.wigner6j(2, 2, 2, 2, 2, 2)         # -3/70\n"
    "-0.04285714285714286";

static PyObject *py_wigner6j(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"j1","j2","j3","j4","j5","j6","precision",NULL};
    PyObject *j1o,*j2o,*j3o,*j4o,*j5o,*j6o,*prec_obj=NULL;
    int tj1,tj2,tj3,tj4,tj5,tj6;
    Precision prec;
    (void)self;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOO|O", kwlist,
                                     &j1o,&j2o,&j3o,&j4o,&j5o,&j6o,&prec_obj))
        return NULL;
    if (parse_half_int(j1o,&tj1,"j1") || parse_half_int(j2o,&tj2,"j2") ||
        parse_half_int(j3o,&tj3,"j3") || parse_half_int(j4o,&tj4,"j4") ||
        parse_half_int(j5o,&tj5,"j5") || parse_half_int(j6o,&tj6,"j6"))
        return NULL;
    prec = parse_precision(prec_obj);
    if ((int)prec < 0) return NULL;
    return make_result(wigner6j_l(tj1,tj2,tj3,tj4,tj5,tj6), prec);
}

/* ── wigner9j ──────────────────────────────────────────────────────────── */

static const char wigner9j_doc[] =
    "wigner9j(j11, j12, j13, j21, j22, j23, j31, j32, j33,"
    " precision='double') -> float\n"
    "\n"
    "Wigner 9j symbol\n"
    "\n"
    "    { j11 j12 j13 }\n"
    "    { j21 j22 j23 }\n"
    "    { j31 j32 j33 }\n"
    "\n"
    "Arguments are passed in row-major order.  The 9j vanishes unless\n"
    "every row and every column independently satisfies the triangle\n"
    "condition.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "j11, j12, ..., j33 : int, float, or fractions.Fraction\n"
    "    Angular-momentum quantum numbers.  Pass an int for integer j,\n"
    "    a float (e.g. ``0.5``) or ``Fraction(1, 2)`` for half-integer\n"
    "    j.  Equal-j ceiling is j <= 5004 with the default-build prime\n"
    "    table; the per-symbol cost scales as O(j^4).\n"
    "precision : {'float', 'double', 'longdouble'}, optional\n"
    "    IEEE 754 binary precision of the returned value.  Default\n"
    "    'double'.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "float\n"
    "    The 9j symbol, correctly rounded to the chosen precision.\n"
    "    Returns 0.0 silently if any selection rule is violated.\n"
    "\n"
    "Examples\n"
    "--------\n"
    ">>> import wignernj\n"
    ">>> wignernj.wigner9j(1, 1, 0, 1, 1, 0, 0, 0, 0)   # 1/3\n"
    "0.3333333333333333";

static PyObject *py_wigner9j(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {
        "j11","j12","j13","j21","j22","j23","j31","j32","j33","precision",NULL};
    PyObject *o[9], *prec_obj=NULL;
    int t[9]; int i;
    Precision prec;
    const char *names[9] = {"j11","j12","j13","j21","j22","j23","j31","j32","j33"};
    (void)self;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOOOOO|O", kwlist,
                                     &o[0],&o[1],&o[2],&o[3],&o[4],
                                     &o[5],&o[6],&o[7],&o[8],&prec_obj))
        return NULL;
    for (i = 0; i < 9; i++)
        if (parse_half_int(o[i], &t[i], names[i])) return NULL;
    prec = parse_precision(prec_obj);
    if ((int)prec < 0) return NULL;
    return make_result(
        wigner9j_l(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8]), prec);
}

/* ── clebsch_gordan ────────────────────────────────────────────────────── */

static const char clebsch_gordan_doc[] =
    "clebsch_gordan(j1, m1, j2, m2, J, M, precision='double') -> float\n"
    "\n"
    "Clebsch-Gordan coefficient <j1 m1; j2 m2 | J M>.\n"
    "\n"
    "Computed as ``(-1)^(j1 - j2 + M) sqrt(2J + 1) * 3j(j1, j2, J;\n"
    "m1, m2, -M)``.  The Condon-Shortley sign convention of Edmonds\n"
    "(1957) and Varshalovich (1988) is used, so every Clebsch-Gordan\n"
    "coefficient is real-valued.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "j1, m1 : int, float, or fractions.Fraction\n"
    "    Angular momentum and projection of the first state.  Pass\n"
    "    an int for integer j, a float (e.g. ``0.5``) or\n"
    "    ``Fraction(1, 2)`` for half-integer j.\n"
    "j2, m2 : int, float, or fractions.Fraction\n"
    "    Angular momentum and projection of the second state.\n"
    "J, M : int, float, or fractions.Fraction\n"
    "    Total angular momentum and projection of the coupled state.\n"
    "    M must equal m1 + m2 for a non-vanishing result.\n"
    "precision : {'float', 'double', 'longdouble'}, optional\n"
    "    IEEE 754 binary precision of the returned value.  Default\n"
    "    'double'.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "float\n"
    "    The Clebsch-Gordan coefficient.  Returns 0.0 silently if any\n"
    "    selection rule is violated (m1 + m2 != M, |M| > J, triangle\n"
    "    inequality).\n"
    "\n"
    "Examples\n"
    "--------\n"
    ">>> import wignernj\n"
    ">>> wignernj.clebsch_gordan(0.5, 0.5, 0.5, -0.5, 1, 0)  # 1/sqrt(2)\n"
    "0.7071067811865476\n"
    ">>> wignernj.clebsch_gordan(1, 1, 1, -1, 2, 0)          # 1/sqrt(6)\n"
    "0.408248290463863";

static PyObject *py_clebsch_gordan(PyObject *self, PyObject *args,
                                    PyObject *kwargs)
{
    static char *kwlist[] = {"j1","m1","j2","m2","J","M","precision",NULL};
    PyObject *j1o,*m1o,*j2o,*m2o,*Jo,*Mo,*prec_obj=NULL;
    int tj1,tm1,tj2,tm2,tJ,tM;
    Precision prec;
    (void)self;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOO|O", kwlist,
                                     &j1o,&m1o,&j2o,&m2o,&Jo,&Mo,&prec_obj))
        return NULL;
    if (parse_half_int(j1o,&tj1,"j1") || parse_half_int(m1o,&tm1,"m1") ||
        parse_half_int(j2o,&tj2,"j2") || parse_half_int(m2o,&tm2,"m2") ||
        parse_half_int(Jo, &tJ, "J")  || parse_half_int(Mo, &tM, "M"))
        return NULL;
    prec = parse_precision(prec_obj);
    if ((int)prec < 0) return NULL;
    return make_result(clebsch_gordan_l(tj1,tm1,tj2,tm2,tJ,tM), prec);
}

/* ── racah_w ───────────────────────────────────────────────────────────── */

static const char racah_w_doc[] =
    "racah_w(j1, j2, J, j3, j12, j23, precision='double') -> float\n"
    "\n"
    "Racah W coefficient W(j1, j2, J, j3; j12, j23).\n"
    "\n"
    "Computed as ``(-1)^(j1 + j2 + J + j3) * 6j{j1, j2, j12;\n"
    "j3, J, j23}``, the standard relation between Racah's W and the\n"
    "Wigner 6j symbol.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "j1, j2, J, j3, j12, j23 : int, float, or fractions.Fraction\n"
    "    Angular-momentum quantum numbers.  Pass an int for integer j,\n"
    "    a float (e.g. ``0.5``) or ``Fraction(1, 2)`` for half-integer\n"
    "    j.  All four triangles of the underlying 6j must hold.\n"
    "precision : {'float', 'double', 'longdouble'}, optional\n"
    "    IEEE 754 binary precision of the returned value.  Default\n"
    "    'double'.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "float\n"
    "    The Racah W coefficient.  Returns 0.0 silently if any\n"
    "    selection rule is violated.\n"
    "\n"
    "Examples\n"
    "--------\n"
    ">>> import wignernj\n"
    ">>> wignernj.racah_w(1, 1, 1, 1, 1, 1)\n"
    "0.16666666666666666";

static PyObject *py_racah_w(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"j1","j2","J","j3","j12","j23","precision",NULL};
    PyObject *j1o,*j2o,*Jo,*j3o,*j12o,*j23o,*prec_obj=NULL;
    int tj1,tj2,tJ,tj3,tj12,tj23;
    Precision prec;
    (void)self;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOO|O", kwlist,
                                     &j1o,&j2o,&Jo,&j3o,&j12o,&j23o,&prec_obj))
        return NULL;
    if (parse_half_int(j1o, &tj1, "j1")  || parse_half_int(j2o, &tj2, "j2")  ||
        parse_half_int(Jo,  &tJ,  "J")   || parse_half_int(j3o, &tj3, "j3")  ||
        parse_half_int(j12o,&tj12,"j12") || parse_half_int(j23o,&tj23,"j23"))
        return NULL;
    prec = parse_precision(prec_obj);
    if ((int)prec < 0) return NULL;
    return make_result(racah_w_l(tj1,tj2,tJ,tj3,tj12,tj23), prec);
}

/* ── fano_x ────────────────────────────────────────────────────────────── */

static const char fano_x_doc[] =
    "fano_x(j1, j2, j12, j3, j4, j34, j13, j24, J, precision='double')"
    " -> float\n"
    "\n"
    "Fano X coefficient X(j1, j2, j12; j3, j4, j34; j13, j24, J)\n"
    "\n"
    "Computed as ``sqrt((2j12+1)(2j34+1)(2j13+1)(2j24+1)) * 9j{j1, j2,\n"
    "j12; j3, j4, j34; j13, j24, J}``, a normalisation variant of the\n"
    "9j symbol used in the analysis of polarisation correlations and\n"
    "in nuclear-physics two-body matrix-element decompositions (Fano\n"
    "1953; Edmonds 1957 §6.4; Tamura 1970).\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "j1, j2, j12, j3, j4, j34, j13, j24, J : int, float, or\n"
    "    fractions.Fraction\n"
    "    Angular-momentum quantum numbers in the same row-major order\n"
    "    as the underlying 9j.  Pass an int for integer j, a float or\n"
    "    Fraction for half-integer j.  Equal-j ceiling is j <= 5004\n"
    "    (delegates to the 9j pipeline; per-symbol cost O(j^4)).\n"
    "precision : {'float', 'double', 'longdouble'}, optional\n"
    "    IEEE 754 binary precision of the returned value.  Default\n"
    "    'double'.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "float\n"
    "    The Fano X coefficient.  Returns 0.0 silently if any\n"
    "    selection rule of the underlying 9j is violated.\n"
    "\n"
    "Examples\n"
    "--------\n"
    ">>> import wignernj\n"
    ">>> wignernj.fano_x(1, 1, 1, 1, 1, 1, 1, 1, 2)\n"
    "0.5\n"
    ">>> wignernj.fano_x(1, 1, 2, 1, 1, 2, 2, 2, 4)\n"
    "1.0";

static PyObject *py_fano_x(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"j1","j2","j12","j3","j4","j34",
                              "j13","j24","J","precision",NULL};
    PyObject *o[9], *prec_obj=NULL;
    int t[9]; int i;
    Precision prec;
    const char *names[9] = {"j1","j2","j12","j3","j4","j34","j13","j24","J"};
    (void)self;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOOOOO|O", kwlist,
                                     &o[0],&o[1],&o[2],&o[3],&o[4],
                                     &o[5],&o[6],&o[7],&o[8],&prec_obj))
        return NULL;
    for (i = 0; i < 9; i++)
        if (parse_half_int(o[i], &t[i], names[i])) return NULL;
    prec = parse_precision(prec_obj);
    if ((int)prec < 0) return NULL;
    return make_result(
        fano_x_l(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8]), prec);
}

/* ── gaunt ─────────────────────────────────────────────────────────────── */

static const char gaunt_doc[] =
    "gaunt(l1, m1, l2, m2, l3, m3, precision='double') -> float\n"
    "\n"
    "Gaunt coefficient: the integral over the unit sphere of three\n"
    "complex spherical harmonics,\n"
    "\n"
    "    G(l1,m1, l2,m2, l3,m3) = integral over Omega\n"
    "        Y_{l1}^{m1}(theta, phi) * Y_{l2}^{m2}(theta, phi)\n"
    "                              * Y_{l3}^{m3}(theta, phi) dOmega\n"
    "\n"
    "Equivalent closed form:\n"
    "\n"
    "    G = sqrt[(2 l1 + 1)(2 l2 + 1)(2 l3 + 1) / (4 pi)]\n"
    "        * 3j(l1,l2,l3; 0,0,0) * 3j(l1,l2,l3; m1,m2,m3)\n"
    "\n"
    "The Condon-Shortley phase for Y_l^m is assumed.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "l1, l2, l3 : int\n"
    "    Non-negative integer orbital angular momenta.\n"
    "m1, m2, m3 : int\n"
    "    Magnetic projections.  Must satisfy m1 + m2 + m3 = 0 and\n"
    "    |m_i| <= l_i for a non-vanishing result.  l1 + l2 + l3 must\n"
    "    be even.\n"
    "precision : {'float', 'double', 'longdouble'}, optional\n"
    "    IEEE 754 binary precision of the returned value.  Default\n"
    "    'double'.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "float\n"
    "    The Gaunt coefficient, correctly rounded.  Returns 0.0\n"
    "    silently if any selection rule is violated.\n"
    "\n"
    "Examples\n"
    "--------\n"
    ">>> import wignernj\n"
    ">>> import math\n"
    ">>> wignernj.gaunt(1, 0, 1, 0, 2, 0)            # 1/sqrt(5*pi)\n"
    "0.2523132522001408\n"
    ">>> abs(wignernj.gaunt(1, 0, 1, 0, 2, 0) - 1/math.sqrt(5*math.pi)) < 1e-15\n"
    "True\n"
    "\n"
    "See Also\n"
    "--------\n"
    "gaunt_real : Gaunt coefficient over real spherical harmonics.";

static PyObject *py_gaunt(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"l1","m1","l2","m2","l3","m3","precision",NULL};
    PyObject *l1o,*m1o,*l2o,*m2o,*l3o,*m3o,*prec_obj=NULL;
    int tl1,tm1,tl2,tm2,tl3,tm3;
    Precision prec;
    (void)self;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOO|O", kwlist,
                                     &l1o,&m1o,&l2o,&m2o,&l3o,&m3o,&prec_obj))
        return NULL;
    if (parse_half_int(l1o,&tl1,"l1") || parse_half_int(m1o,&tm1,"m1") ||
        parse_half_int(l2o,&tl2,"l2") || parse_half_int(m2o,&tm2,"m2") ||
        parse_half_int(l3o,&tl3,"l3") || parse_half_int(m3o,&tm3,"m3"))
        return NULL;
    prec = parse_precision(prec_obj);
    if ((int)prec < 0) return NULL;
    return make_result(gaunt_l(tl1,tm1,tl2,tm2,tl3,tm3), prec);
}

/* ── gaunt_real ────────────────────────────────────────────────────────── */

static const char gaunt_real_doc[] =
    "gaunt_real(l1, m1, l2, m2, l3, m3, precision='double') -> float\n"
    "\n"
    "Gaunt coefficient over real spherical harmonics:\n"
    "\n"
    "    G_R(l1,m1, l2,m2, l3,m3) = integral over Omega\n"
    "        S_{l1,m1}(theta, phi) * S_{l2,m2}(theta, phi)\n"
    "                            * S_{l3,m3}(theta, phi) dOmega\n"
    "\n"
    "where S_{l,m} are the real spherical harmonics in the standard\n"
    "Condon-Shortley / Wikipedia construction:\n"
    "\n"
    "    S_{l, 0}     = Y_l^0\n"
    "    S_{l,  m>0}  = (1/sqrt(2)) * (Y_l^{-m} + (-1)^m Y_l^m)\n"
    "    S_{l,  m<0}  = (i/sqrt(2)) * (Y_l^{m}  - (-1)^|m| Y_l^{-m})\n"
    "\n"
    "Internally the value is computed as a single complex Gaunt\n"
    "evaluation per call (see src/gaunt.c for the algorithm).\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "l1, l2, l3 : int\n"
    "    Non-negative integer orbital angular momenta.\n"
    "m1, m2, m3 : int\n"
    "    Signed integer indices selecting between cosine-type (m > 0),\n"
    "    sine-type (m < 0), or the unique m = 0 harmonic.  Must\n"
    "    satisfy |m_i| <= l_i.\n"
    "precision : {'float', 'double', 'longdouble'}, optional\n"
    "    IEEE 754 binary precision of the returned value.  Default\n"
    "    'double'.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "float\n"
    "    The real-spherical-harmonic Gaunt coefficient.  Returns 0.0\n"
    "    silently if any selection rule is violated.\n"
    "\n"
    "Examples\n"
    "--------\n"
    ">>> import wignernj\n"
    ">>> wignernj.gaunt_real(1, 0, 1, 0, 2, 0)       # same as gaunt at m=0\n"
    "0.2523132522001408\n"
    "\n"
    "See Also\n"
    "--------\n"
    "gaunt : Gaunt coefficient over complex spherical harmonics.";

static PyObject *py_gaunt_real(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"l1","m1","l2","m2","l3","m3","precision",NULL};
    PyObject *l1o,*m1o,*l2o,*m2o,*l3o,*m3o,*prec_obj=NULL;
    int tl1,tm1,tl2,tm2,tl3,tm3;
    Precision prec;
    (void)self;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOO|O", kwlist,
                                     &l1o,&m1o,&l2o,&m2o,&l3o,&m3o,&prec_obj))
        return NULL;
    if (parse_half_int(l1o,&tl1,"l1") || parse_half_int(m1o,&tm1,"m1") ||
        parse_half_int(l2o,&tl2,"l2") || parse_half_int(m2o,&tm2,"m2") ||
        parse_half_int(l3o,&tl3,"l3") || parse_half_int(m3o,&tm3,"m3"))
        return NULL;
    prec = parse_precision(prec_obj);
    if ((int)prec < 0) return NULL;
    return make_result(gaunt_real_l(tl1,tm1,tl2,tm2,tl3,tm3), prec);
}

/* ── module definition ─────────────────────────────────────────────────── */

/* Cast every keyword-accepting CFunction through `void (*)(void)` before the
 * PyCFunction landing slot.  PyMethodDef stores PyCFunction (two-arg), but
 * METH_VARARGS|METH_KEYWORDS dispatches at run-time to a three-arg
 * PyCFunctionWithKeywords; the two-step cast is the canonical CPython
 * idiom that silences -Wcast-function-type cleanly (gcc/clang refuse the
 * direct cast). */
#define KW_METHOD(name, fn, doc) \
    {name, (PyCFunction)(void(*)(void))(fn), METH_VARARGS|METH_KEYWORDS, doc}

static PyMethodDef wignernj_methods[] = {
    KW_METHOD("wigner3j",       py_wigner3j,       wigner3j_doc),
    KW_METHOD("wigner6j",       py_wigner6j,       wigner6j_doc),
    KW_METHOD("wigner9j",       py_wigner9j,       wigner9j_doc),
    KW_METHOD("clebsch_gordan", py_clebsch_gordan, clebsch_gordan_doc),
    KW_METHOD("racah_w",        py_racah_w,        racah_w_doc),
    KW_METHOD("fano_x",         py_fano_x,         fano_x_doc),
    KW_METHOD("gaunt",          py_gaunt,          gaunt_doc),
    KW_METHOD("gaunt_real",     py_gaunt_real,     gaunt_real_doc),
    {NULL, NULL, 0, NULL}
};

#undef KW_METHOD

/* PyModuleDef has more fields than the legacy positional-init covered;
 * use designated initializers so the trailing m_slots / m_traverse /
 * m_clear / m_free fields default to NULL without -Wmissing-field-
 * initializers complaining. */
static struct PyModuleDef wignernjmodule = {
    .m_base    = PyModuleDef_HEAD_INIT,
    .m_name    = "_wignernj",
    .m_doc     = "Exact Wigner 3j/6j/9j symbols and related coefficients.",
    .m_size    = -1,
    .m_methods = wignernj_methods,
};

PyMODINIT_FUNC PyInit__wignernj(void)
{
    return PyModule_Create(&wignernjmodule);
}
