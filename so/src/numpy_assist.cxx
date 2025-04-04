#define NO_IMPORT_ARRAY

#include <so/numpy_assist.h>


template <>
bool check_buffer_type<float>(const Py_buffer &view) {
    return _check_buffer_helper<float>(view, "f");
}

template <>
bool check_buffer_type<double>(const Py_buffer &view) {
    return _check_buffer_helper<double>(view, "d");
}

template <>
std::string type_name<int64_t>() {
    return "int64";
}

template <>
std::string type_name<int32_t>() {
    return "int32";
}

template <>
std::string type_name<float>() {
    return "float32";
}

template <>
std::string type_name<double>() {
    return "float64";
}

// Work around unused function warning.  Define a function that is never called.
void dummy_function_use() {
    std::string ret;
    ret = type_name<int32_t>();
    ret = type_name<int64_t>();
    ret = type_name<float>();
    ret = type_name<double>();

    auto p = (Py_buffer*)calloc(1, sizeof(Py_buffer));
    auto view = std::shared_ptr<Py_buffer>(p, PyBuffer_Release);
    bool check;
    check = check_buffer_type<float>(*view.get());
    check = check_buffer_type<double>(*view.get());
    (void)check;
    return;
}

