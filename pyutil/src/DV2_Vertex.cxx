#include <Python.h>

#include "WCPPyUtil/DV2_Vertex.h"

#include <cstring>
#include <stdexcept>
#include <vector>

// Same embedding pattern as SCN_Vertex.cxx: initialise the interpreter (safe
// to repeat), import the module by name, call one function with bytes
// buffers, read a bytes result.  Py_Finalize() is never called (torch).

namespace {
    std::string fetch_python_error()
    {
        std::string pyerr;
        if (PyErr_Occurred()) {
            PyObject *ptype, *pvalue, *ptraceback;
            PyErr_Fetch(&ptype, &pvalue, &ptraceback);
            PyErr_NormalizeException(&ptype, &pvalue, &ptraceback);
            PyObject* pstr = pvalue ? PyObject_Str(pvalue) : nullptr;
            if (pstr) {
                const char* s = PyUnicode_AsUTF8(pstr);
                if (s) pyerr = std::string(": ") + s;
                Py_DECREF(pstr);
            }
            Py_XDECREF(ptype);
            Py_XDECREF(pvalue);
            Py_XDECREF(ptraceback);
        }
        return pyerr;
    }

    PyObject* bytes_of(const std::vector<WCPPyUtil::FLOAT>& v)
    {
        PyObject* b = PyBytes_FromStringAndSize(reinterpret_cast<const char*>(v.data()),
                                                v.size() * sizeof(WCPPyUtil::FLOAT));
        if (!b || !PyBytes_Check(b)) {
            Py_XDECREF(b);
            throw std::runtime_error("DV2_Vertex: failed to create a bytes buffer");
        }
        return b;
    }
}

namespace WCPPyUtil {
    std::vector<FLOAT> DV2_Vertex(const std::string& module, const std::string& function,
                                  const std::string& weights,
                                  const std::vector<std::vector<FLOAT> >& input, int n_vertex_rows,
                                  const std::vector<FLOAT>& trad, double q_scale, double q_offset,
                                  const std::string& dtype)
    {
        if (input.size() != 4) throw std::runtime_error("DV2_Vertex: input.size() != 4");
        const size_t npts = input[3].size();
        for (size_t i = 0; i < 3; ++i) {
            if (input[i].size() != npts) throw std::runtime_error("DV2_Vertex: input size mismatch");
        }
        if (trad.size() != 3) throw std::runtime_error("DV2_Vertex: trad must have 3 values");
        if (n_vertex_rows < 0 || static_cast<size_t>(n_vertex_rows) > npts) {
            throw std::runtime_error("DV2_Vertex: n_vertex_rows out of range");
        }

        Py_Initialize();

        PyObject* pName = PyUnicode_FromString(module.c_str());
        PyObject* pModule = PyImport_Import(pName);
        Py_DECREF(pName);
        if (!pModule) {
            throw std::runtime_error("DV2_Vertex: import failed for module: " + module + fetch_python_error());
        }
        PyObject* pFunc = PyObject_GetAttrString(pModule, function.c_str());   // new reference
        if (!pFunc || !PyCallable_Check(pFunc)) {
            Py_XDECREF(pFunc);
            Py_DECREF(pModule);
            throw std::runtime_error("DV2_Vertex: function not found or not callable: " + function);
        }

        PyObject* args[10] = {nullptr};
        args[0] = PyUnicode_FromString(weights.c_str());
        args[1] = bytes_of(input[0]);
        args[2] = bytes_of(input[1]);
        args[3] = bytes_of(input[2]);
        args[4] = bytes_of(input[3]);
        args[5] = PyLong_FromLong(n_vertex_rows);
        args[6] = bytes_of(trad);
        args[7] = PyFloat_FromDouble(q_scale);
        args[8] = PyFloat_FromDouble(q_offset);
        args[9] = PyUnicode_FromString(dtype.c_str());

        PyObject* pValue = PyObject_CallFunctionObjArgs(pFunc, args[0], args[1], args[2], args[3], args[4],
                                                        args[5], args[6], args[7], args[8], args[9], NULL);
        for (auto* a : args) Py_XDECREF(a);
        Py_DECREF(pFunc);
        Py_DECREF(pModule);

        if (!pValue) {
            throw std::runtime_error("DV2_Vertex: Python function call failed" + fetch_python_error());
        }
        if (!PyBytes_Check(pValue)) {
            Py_DECREF(pValue);
            throw std::runtime_error("DV2_Vertex: Python function did not return bytes");
        }
        const size_t ret_size = PyBytes_Size(pValue);
        if (ret_size % sizeof(FLOAT) != 0) {
            Py_DECREF(pValue);
            throw std::runtime_error("DV2_Vertex: returned bytes are not a whole number of floats");
        }
        std::vector<FLOAT> ret(ret_size / sizeof(FLOAT));
        if (ret_size) std::memcpy(ret.data(), PyBytes_AsString(pValue), ret_size);
        Py_DECREF(pValue);
        return ret;
    }
}  // namespace WCPPyUtil
