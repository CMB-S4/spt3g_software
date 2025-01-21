#pragma once

#include <boost/python.hpp>
#include <exception>

// so3g_exception is our internal base class, which defines the
// interface we use for converting C++ exceptions to python.

class so3g_exception : std::exception
{
public:
    std::string text;
    so3g_exception() {};

    so3g_exception(std::string text) :
        text{text} {}

    virtual std::string msg_for_python() const throw() {
        return text;
    }
};


// The base classes here are mapped to specific Python exceptions.
// They are registered with boost python in exceptions.cxx.  Sure, you
// can use these directly.  Why not.

class RuntimeError_exception : public so3g_exception {
    using so3g_exception::so3g_exception;
};
class TypeError_exception : public so3g_exception {
    using so3g_exception::so3g_exception;
};
class ValueError_exception : public so3g_exception {
    using so3g_exception::so3g_exception;
};


// The exceptions below should be used when processing objects with
// the buffer protocol (probably numpy arrays).

class buffer_exception : public TypeError_exception
{
public:
    std::string var_name;
    buffer_exception(std::string var_name) : var_name{var_name} {}

    std::string msg_for_python() const throw() {
        std::ostringstream s;
        s << "Argument '" << var_name << "' does not expose buffer protocol, "
            "is not contiguous, or does not export a format.";
        return s.str();
    }
};

class shape_exception : public RuntimeError_exception
{
public:
    std::string var_name;
    std::string detail;
    shape_exception(std::string var_name, std::string detail) :
        var_name{var_name}, detail(detail) {}

    std::string msg_for_python() const throw() {
        std::ostringstream s;
        s << "Buffer '" << var_name << "' has incompatible shape: "
          << detail << ".";
        return s.str();
    }
};

class dtype_exception : public ValueError_exception
{
public:
    std::string var_name;
    std::string type_str;
    dtype_exception(std::string var_name, std::string type_str) :
        var_name{var_name}, type_str{type_str} {}

    std::string msg_for_python() const throw() {
        std::ostringstream s;
        s << "Expected buffer '" << var_name << "' to contain items of type "
          << type_str << ".";
        return s.str();
    }
};

class agreement_exception : public RuntimeError_exception
{
public:
    std::string var1, var2, prop;
    agreement_exception(std::string var1, std::string var2, std::string prop) :
        var1{var1}, var2{var2}, prop{prop} {}

    std::string msg_for_python() const throw() {
        std::ostringstream s;
        s << "Expected buffers '" << var1 << "' and '" << var2 << "' to have "
          << "the same " << prop << ".";
        return s.str();
    }
};

class tiling_exception : public RuntimeError_exception
{
public:
    int tile_idx;
    std::string msg;
    tiling_exception(int tile_idx, std::string msg) :
        tile_idx{tile_idx}, msg{msg} {}

    std::string msg_for_python() const throw() {
        std::ostringstream s;
        s << "Tiling problem (index " << tile_idx << "): " << msg;
        return s.str();
    }
};

class general_agreement_exception : public ValueError_exception
{
public:
    std::string text;
    general_agreement_exception(std::string text) :
        text{text} {}

    std::string msg_for_python() const throw() {
        return text;
    }
};
