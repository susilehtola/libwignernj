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
#include "wigner.h"

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
    "wigner3j(j1, j2, j3, m1, m2, m3, precision='double') -> float\n\n"
    "Wigner 3j symbol  ( j1  j2  j3 )\n"
    "                  ( m1  m2  m3 )\n\n"
    "Arguments may be int, float (half-integer), or fractions.Fraction.\n"
    "Returns 0.0 if selection rules are violated (not an error).\n"
    "precision: 'float', 'double' (default), or 'longdouble'.";

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
    "wigner6j(j1, j2, j3, j4, j5, j6, precision='double') -> float\n\n"
    "Wigner 6j symbol  { j1 j2 j3 }\n"
    "                  { j4 j5 j6 }";

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
    "wigner9j(j11,j12,j13, j21,j22,j23, j31,j32,j33, precision='double')"
    " -> float\n\n"
    "Wigner 9j symbol (row-major order).";

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
    "clebsch_gordan(j1, m1, j2, m2, J, M, precision='double') -> float\n\n"
    "Clebsch-Gordan coefficient <j1 m1; j2 m2 | J M>.";

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
    "racah_w(j1, j2, J, j3, j12, j23, precision='double') -> float\n\n"
    "Racah W-coefficient W(j1 j2 J j3; j12 j23).";

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
    " -> float\n\n"
    "Fano X-coefficient X(j1 j2 j12; j3 j4 j34; j13 j24 J).\n"
    "Equals sqrt[(2j12+1)(2j34+1)(2j13+1)(2j24+1)] times the corresponding\n"
    "9j symbol.";

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
    "gaunt(l1, m1, l2, m2, l3, m3, precision='double') -> float\n\n"
    "Gaunt coefficient: integral Y_{l1}^{m1} Y_{l2}^{m2} Y_{l3}^{m3} dΩ.\n"
    "l arguments must be non-negative integers; m1+m2+m3 must equal 0.";

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
    "gaunt_real(l1, m1, l2, m2, l3, m3, precision='double') -> float\n\n"
    "Real-spherical-harmonic Gaunt coefficient:\n"
    "    integral S_{l1,m1} S_{l2,m2} S_{l3,m3} dΩ\n"
    "with the standard Condon-Shortley / Wikipedia convention for the real\n"
    "spherical harmonics S_{l,m}.  l must be non-negative integers and\n"
    "|m_i| <= l_i; signed m selects cosine-type (m>0), sine-type (m<0), or\n"
    "the m=0 harmonic.";

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

static PyMethodDef wigner_methods[] = {
    {"wigner3j",       (PyCFunction)py_wigner3j,       METH_VARARGS|METH_KEYWORDS, wigner3j_doc},
    {"wigner6j",       (PyCFunction)py_wigner6j,       METH_VARARGS|METH_KEYWORDS, wigner6j_doc},
    {"wigner9j",       (PyCFunction)py_wigner9j,       METH_VARARGS|METH_KEYWORDS, wigner9j_doc},
    {"clebsch_gordan", (PyCFunction)py_clebsch_gordan, METH_VARARGS|METH_KEYWORDS, clebsch_gordan_doc},
    {"racah_w",        (PyCFunction)py_racah_w,        METH_VARARGS|METH_KEYWORDS, racah_w_doc},
    {"fano_x",         (PyCFunction)py_fano_x,         METH_VARARGS|METH_KEYWORDS, fano_x_doc},
    {"gaunt",          (PyCFunction)py_gaunt,          METH_VARARGS|METH_KEYWORDS, gaunt_doc},
    {"gaunt_real",     (PyCFunction)py_gaunt_real,     METH_VARARGS|METH_KEYWORDS, gaunt_real_doc},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef wignermodule = {
    PyModuleDef_HEAD_INIT,
    "_wigner",
    "Exact Wigner 3j/6j/9j symbols and related coefficients.",
    -1,
    wigner_methods
};

PyMODINIT_FUNC PyInit__wigner(void)
{
    return PyModule_Create(&wignermodule);
}
