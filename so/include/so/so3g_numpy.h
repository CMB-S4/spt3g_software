#pragma once

// Source files, except for main.cxx, need to define NO_IMPORT_ARRAY
// before importing this file.
//
// Some macro defines are necessary when working with numpy C api
// across multiple translation units (i.e. source files):
//
// - https://sourceforge.net/p/numpy/mailman/message/5700519/
// - https://stackoverflow.com/questions/38003707/trouble-with-numpy-c-api
//
// As a result, independently of the boost numpy stuff, we need to
// tell numpy to drop in a single copy of the initialized API function
// pointers.  All files using the API should define the
// PY_ARRAY_UNIQUE_SYMBOL variable to the same value (by including
// this header file).  Then we call import_array() in one of those
// source files (main.cxx), and define NO_IMPORT_ARRAY in the others.
// Then, from this header file, it is safe to include arrayobject.h.

#define PY_ARRAY_UNIQUE_SYMBOL Py_Array_API_SO3G
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

