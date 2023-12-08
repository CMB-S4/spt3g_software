#include <dlfcn.h>

#include <Python.h>

static PyObject* load_dynamic_library(PyObject* mod, PyObject* args, PyObject* kwds){
	static const char* kwlist[] = {"full_name", "name", "path", NULL};
	char* full_name=NULL;
	char* name=NULL;
	char* path=NULL;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "sss", (char**)kwlist, &full_name, &name, &path))
		return Py_None;
	
	void* h = dlopen(path, RTLD_NOW | RTLD_GLOBAL);
    char* errmsg = dlerror();
    
    if(h == NULL || errmsg != NULL){
    	char err_buf[1024];
    	snprintf(err_buf, 1024, "dynamic loading error: loading '%s' from '%s': %s", name, path, 
    	         (errmsg == NULL ? "dlopen: unknown error" : errmsg));
    	PyErr_SetString(PyExc_RuntimeError, err_buf);
		return Py_None;
    }
    PyObject* (*init)(void);
    char initName[256];
    int print_res = 0;
#if PY_MAJOR_VERSION >= 3
	print_res = snprintf(initName, 256, "PyInit_%s", name);
#else
	print_res = snprintf(initName, 256, "init%s", name);
#endif
	if(print_res<0){
		char err_buf[1024];
    	snprintf(err_buf, 1024, "dynamic loading error: loading '%s' from '%s': "
    	         "module init function name too long", name, path);
    	PyErr_SetString(PyExc_RuntimeError, err_buf);
		return Py_None;
	}

    *(void **) (&init) = dlsym(h,initName);
    errmsg = dlerror();
    if(errmsg != NULL){
    	char err_buf[1024];
    	snprintf(err_buf, 1024, "dynamic loading error: loading '%s' from '%s': %s", name, path,
    	         (errmsg == NULL ? "dlsym: unknown error" : errmsg));
    	printf("%s\n", err_buf);
    	PyErr_SetString(PyExc_RuntimeError, err_buf);
		return Py_None;
    }
    
#if PY_MAJOR_VERSION >= 3
    PyObject* module = (*init)();
    PyObject_SetAttrString(module, "__file__", PyUnicode_FromString(path));
#else
	(*init)();
	PyObject* module = PyDict_GetItemString(PyImport_GetModuleDict(), full_name);
	if(nmod==NULL){
		char err_buf[1024];
		snprintf(err_buf, 1024, "dynamic loading error: %s not in global module dict", full_name);
		PyErr_SetString(PyExc_RuntimeError, err_buf);
		return Py_None;
	}
	PyObject_SetAttrString(module, "__file__", PyString_FromString(path));
#endif
	return module;
}

static PyMethodDef spt3g_dload_methods[] = {
	{"load_dynamic", (PyCFunction)&load_dynamic_library, METH_VARARGS | METH_KEYWORDS, 
	 "\n"
	 "Load a compiled extension module. No restriction is placed on the relationship\n"
	 "between the name of the module and the name of the dynamic library file which\n"
	 "implements it; e.g. module `foo` need not reside in `foo.so`.\n"
	 "\n"
	 "Arguments\n"
	 "---------\n"
	 "full_name : str\n"
	 "    The module's fully qualified name\n"
	 "name : str\n"
	 "    The module's name\n"
	 "path : str\n"
	 "    The path to the dynamic library which implements the module\n"
	 "\n"
	 "Returns\n"
	 "-------\n"
	 "The loaded module\n"
	},
	{NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"spt3g_dload",
	"A beach-head module for loading other compiled spt3g modules",
	0,
	spt3g_dload_methods,
	NULL,
	NULL,
	NULL,
	NULL,
};
#endif

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
PyInit_dload(void){
#else
initdload(void){
#endif
	PyObject* module;
	
#if PY_MAJOR_VERSION >= 3
	module = PyModule_Create(&moduledef);
	PyObject_SetAttrString(module, "__version__", PyUnicode_FromString("1.0.0"));
	return module;
#else
	module = Py_InitModule3("spt3g_dload", spt3g_dload_methods,
	                        "A beach-head module for loading other compiled spt3g modules");
	PyObject_SetAttrString(module, "__version__", PyString_FromString("1.0.0"));
	PyDict_SetItemString(PyImport_GetModuleDict(),_Py_PackageContext,module);
#endif
}