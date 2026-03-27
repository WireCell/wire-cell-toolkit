#include <Python.h>

#include "WCPPyUtil/SCN_Vertex.h"

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <cstring>

#define DebugVar(x) std::cout << "DebugVar: " << #x << ": " << x << std::endl
#define DebugInf(x) std::cout << "DebugInf: " << x << std::endl

template <class VType>
static void print_vec(const std::vector<VType> &data)
{
    for (auto v : data) std::cout << v << " ";
    std::cout << std::endl;
}

namespace WCPPyUtil {
    std::vector<FLOAT> SCN_Vertex(const std::string &module, const std::string &function, const std::string &weights,
                                  const std::vector< std::vector<FLOAT> > &input, const std::string &dtype, const bool verbose)
    {
        if (input.size() != 4) {
            throw std::runtime_error("input.size() != 4");
        }
        auto x = input[0];
        auto y = input[1];
        auto z = input[2];
        auto q = input[3];
        if (verbose) {
            DebugInf("SCN_Vertex: start");
            DebugVar(module);
            DebugVar(function);
            DebugVar(weights);
            DebugVar(dtype);
        }

        size_t npts = q.size();
        if (x.size() != npts || y.size() != npts || z.size() != npts) {
            throw std::runtime_error("input size unmatch");
        }

        PyObject *pName, *pWeights, *pDtype, *pModule, *pDict, *pFunc, *pValue;

        // Initialize the Python interpreter (safe to call multiple times)
        Py_Initialize();

        pName    = PyUnicode_FromString(module.c_str());
        pWeights = PyUnicode_FromString(weights.c_str());
        pDtype   = PyUnicode_FromString(dtype.c_str());

        // Load the module object
        pModule = PyImport_Import(pName);
        if (!pModule) {
            PyErr_Print();
            throw std::runtime_error("SCN_Vertex: import failed for module: " + module);
        }

        // pDict is a borrowed reference
        pDict = PyModule_GetDict(pModule);
        if (!pDict) {
            throw std::runtime_error("SCN_Vertex: PyModule_GetDict failed");
        }

        // pFunc is also a borrowed reference
        pFunc = PyDict_GetItemString(pDict, function.c_str());
        if (!pFunc) {
            throw std::runtime_error("SCN_Vertex: function not found: " + function);
        }

        if (PyCallable_Check(pFunc) != true) {
            throw std::runtime_error("SCN_Vertex: function not callable: " + function);
        }

        size_t input_size = npts * sizeof(FLOAT);
        auto pX = PyBytes_FromStringAndSize((const char *) x.data(), input_size);
        auto pY = PyBytes_FromStringAndSize((const char *) y.data(), input_size);
        auto pZ = PyBytes_FromStringAndSize((const char *) z.data(), input_size);
        auto pQ = PyBytes_FromStringAndSize((const char *) q.data(), input_size);
        if (!PyBytes_Check(pX) || !PyBytes_Check(pY) || !PyBytes_Check(pZ) || !PyBytes_Check(pQ)) {
            throw std::runtime_error("SCN_Vertex: failed to create byte buffers");
        }

        pValue = PyObject_CallFunctionObjArgs(pFunc, pWeights, pX, pY, pZ, pQ, pDtype, NULL);
        if (verbose) {
            DebugInf("PyObject_CallFunctionObjArgs: OK");
        }

        Py_DECREF(pX);
        Py_DECREF(pY);
        Py_DECREF(pZ);
        Py_DECREF(pQ);

        if (!pValue) {
            PyErr_Print();
            throw std::runtime_error("SCN_Vertex: Python function call failed");
        }

        size_t ret_size = PyBytes_Size(pValue);
        if (verbose) {
            DebugVar(ret_size);
        }
        assert(ret_size % sizeof(FLOAT) == 0);
        std::vector<FLOAT> ret(ret_size / sizeof(FLOAT));
        memcpy((char *) ret.data(), (char *) PyBytes_AsString(pValue), ret_size);

        if (verbose) {
            printf("SCN_Vertex: Return of call: ");
            print_vec<FLOAT>(ret);
        }
        Py_DECREF(pValue);

        // Clean up
        Py_DECREF(pModule);
        Py_DECREF(pName);
        Py_DECREF(pWeights);
        Py_DECREF(pDtype);

        // Note: Py_Finalize() is intentionally NOT called here.
        // Calling it causes a segfault when sparseconvnet/torch are loaded.
        // The interpreter stays alive for the process lifetime.

        return ret;
    }
}  // namespace WCPPyUtil
