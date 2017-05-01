#pragma once
#include <pybindings.h>

void coord_gal_to_eq(double *l_in, double *b_in, double *ra_out, double *dec_out,
                     size_t size);

template <typename T>
boost::python::tuple coord_gal_to_eq_py(T l_in, T b_in);
