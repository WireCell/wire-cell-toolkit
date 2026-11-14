#include <Python.h>

#include "WCPPyUtil/DV2_Vertex.h"

#include "GIL.h"

#include <cstring>
#include <stdexcept>
#include <vector>

// Same embedding pattern as SCN_Vertex.cxx: start the interpreter once and take
// the GIL for the call (GIL.h), import the module by name, call one function with bytes
// buffers, read a bytes result.  Py_Finalize() is never called (torch).

namespace {
    std::string fetch_python_error()
    {
        std::string pyerr;
        if (PyErr_Occurred()) {
            PyObject *ptype, *pvalue, *ptraceback;
            PyErr_Fetch(&ptype, &pvalue, &ptraceback);
            PyErr_NormalizeException(&ptype, &pvalue, &ptraceback);
            PyObject* pname = ptype ? PyObject_GetAttrString(ptype, "__name__") : nullptr;
            const char* n = pname ? PyUnicode_AsUTF8(pname) : nullptr;
            PyObject* pstr = pvalue ? PyObject_Str(pvalue) : nullptr;
            const char* s = pstr ? PyUnicode_AsUTF8(pstr) : nullptr;
            pyerr = std::string(": ") + (n ? n : "Python error") + (s && *s ? std::string(": ") + s : std::string());
            Py_XDECREF(pstr);
            Py_XDECREF(pname);
            PyErr_Clear();
            Py_XDECREF(ptype);
            Py_XDECREF(pvalue);
            Py_XDECREF(ptraceback);
        }
        return pyerr;
    }

    // New reference to a bytes copy of v (nullptr if the allocation fails).
    PyObject* bytes_of(const std::vector<WCPPyUtil::FLOAT>& v)
    {
        return PyBytes_FromStringAndSize(reinterpret_cast<const char*>(v.data()),
                                         v.size() * sizeof(WCPPyUtil::FLOAT));
    }

    // The cluster ids travel as int32 (DV2_Vertex.py reads them with np.int32).
    static_assert(sizeof(int) == 4, "DV2_Vertex: space_cluster_id is passed as 4-byte ints");
    PyObject* bytes_of(const std::vector<int>& v)
    {
        return PyBytes_FromStringAndSize(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(int));
    }

    // Owns new references: Py_XDECREF of each on scope exit, also when an
    // exception unwinds, so no error path leaks.  Capacity is reserved up
    // front, so own() never allocates (and never throws) for <= kMax refs.
    class PyRefs {
      public:
        static constexpr size_t kMax = 24;
        PyRefs() { m_refs.reserve(kMax); }
        PyRefs(const PyRefs&) = delete;
        PyRefs& operator=(const PyRefs&) = delete;
        ~PyRefs()
        {
            for (auto it = m_refs.rbegin(); it != m_refs.rend(); ++it) Py_XDECREF(*it);
        }
        PyObject* own(PyObject* o)
        {
            if (m_refs.size() == kMax) {
                Py_XDECREF(o);
                throw std::logic_error("PyRefs: more than kMax references");
            }
            m_refs.push_back(o);
            return o;
        }

      private:
        std::vector<PyObject*> m_refs;
    };

    // The module `module` imported (owned by refs); throws if the import fails.
    PyObject* import_module(PyRefs& refs, const std::string& module, const std::string& who)
    {
        PyObject* pName = refs.own(PyUnicode_FromString(module.c_str()));
        PyObject* pModule = pName ? refs.own(PyImport_Import(pName)) : nullptr;
        if (!pModule) {
            throw std::runtime_error(who + ": import failed for module " + module + fetch_python_error());
        }
        return pModule;
    }

    // module.function (owned by refs); throws if it is missing or not callable.
    PyObject* get_function(PyRefs& refs, PyObject* pModule, const std::string& module, const std::string& function,
                           const std::string& who)
    {
        PyObject* pFunc = refs.own(PyObject_GetAttrString(pModule, function.c_str()));
        if (!pFunc || !PyCallable_Check(pFunc)) {
            PyErr_Clear();
            throw std::runtime_error(who + ": function not found or not callable: " + module + "." + function);
        }
        return pFunc;
    }
}

namespace WCPPyUtil {
    std::string DV2_Load(const std::string& module, const std::string& function, const std::string& weights)
    {
        GIL gil;       // the interpreter is started once; the GIL is held for this call only
        PyRefs refs;   // declared after gil: released while the GIL is still held
        PyObject* pModule = import_module(refs, module, "DV2_Load");
        PyObject* pFunc = get_function(refs, pModule, module, function, "DV2_Load");
        PyObject* pArg = refs.own(PyUnicode_FromString(weights.c_str()));
        if (!pArg) {
            throw std::runtime_error("DV2_Load: could not build the Python argument" + fetch_python_error());
        }
        PyObject* pValue = refs.own(PyObject_CallFunctionObjArgs(pFunc, pArg, NULL));
        if (!pValue) {
            throw std::runtime_error("DV2_Load: loading " + weights + " failed" + fetch_python_error());
        }
        std::string ret;
        if (PyBytes_Check(pValue)) ret.assign(PyBytes_AsString(pValue), PyBytes_Size(pValue));
        return ret;
    }

    // Both DV2_Vertex overloads: space / space_cluster_id are nullptr for the cloud-only call.
    static std::vector<FLOAT> dv2_call(const std::string& module, const std::string& function,
                                       const std::string& weights,
                                       const std::vector<std::vector<FLOAT> >& input, int n_vertex_rows,
                                       const std::vector<FLOAT>& trad, double q_scale, double q_offset,
                                       const std::string& dtype,
                                       const std::vector<std::vector<FLOAT> >* space,
                                       const std::vector<int>* space_cluster_id)
    {
        if (input.size() != 4) throw std::runtime_error("DV2_Vertex: input.size() != 4");
        if (space) {
            if (space->size() != 4) throw std::runtime_error("DV2_Vertex: space.size() != 4");
            const size_t ns = (*space)[3].size();
            for (size_t i = 0; i < 3; ++i) {
                if ((*space)[i].size() != ns) throw std::runtime_error("DV2_Vertex: space size mismatch");
            }
            if (!space_cluster_id || space_cluster_id->size() != ns) {
                throw std::runtime_error("DV2_Vertex: space_cluster_id size mismatch");
            }
        }
        const size_t npts = input[3].size();
        for (size_t i = 0; i < 3; ++i) {
            if (input[i].size() != npts) throw std::runtime_error("DV2_Vertex: input size mismatch");
        }
        if (trad.size() != 3) throw std::runtime_error("DV2_Vertex: trad must have 3 values");
        if (n_vertex_rows < 0 || static_cast<size_t>(n_vertex_rows) > npts) {
            throw std::runtime_error("DV2_Vertex: n_vertex_rows out of range");
        }

        GIL gil;       // the interpreter is started once; the GIL is held for this call only
        PyRefs refs;   // every reference below is released on return or throw, with the GIL held
        PyObject* pModule = import_module(refs, module, "DV2_Vertex");
        PyObject* pFunc = get_function(refs, pModule, module, function, "DV2_Vertex");

        // braced initialisation runs left to right
        PyObject* args[10] = {refs.own(PyUnicode_FromString(weights.c_str())),
                              refs.own(bytes_of(input[0])),
                              refs.own(bytes_of(input[1])),
                              refs.own(bytes_of(input[2])),
                              refs.own(bytes_of(input[3])),
                              refs.own(PyLong_FromLong(n_vertex_rows)),
                              refs.own(bytes_of(trad)),
                              refs.own(PyFloat_FromDouble(q_scale)),
                              refs.own(PyFloat_FromDouble(q_offset)),
                              refs.own(PyUnicode_FromString(dtype.c_str()))};
        // A NULL argument would end PyObject_CallFunctionObjArgs' list early.
        for (auto* a : args) {
            if (!a) throw std::runtime_error("DV2_Vertex: could not build the Python arguments" + fetch_python_error());
        }

        PyObject* pValue = nullptr;
        if (space) {
            PyObject* sargs[5] = {refs.own(bytes_of((*space)[0])),
                                  refs.own(bytes_of((*space)[1])),
                                  refs.own(bytes_of((*space)[2])),
                                  refs.own(bytes_of((*space)[3])),
                                  refs.own(bytes_of(*space_cluster_id))};
            for (auto* a : sargs) {
                if (!a) throw std::runtime_error("DV2_Vertex: could not build the Python arguments" + fetch_python_error());
            }
            pValue = refs.own(PyObject_CallFunctionObjArgs(pFunc, args[0], args[1], args[2], args[3], args[4],
                                                           args[5], args[6], args[7], args[8], args[9],
                                                           sargs[0], sargs[1], sargs[2], sargs[3], sargs[4], NULL));
        }
        else {
            pValue = refs.own(PyObject_CallFunctionObjArgs(pFunc, args[0], args[1], args[2], args[3], args[4],
                                                           args[5], args[6], args[7], args[8], args[9], NULL));
        }
        if (!pValue) {
            throw std::runtime_error("DV2_Vertex: Python function call failed" + fetch_python_error());
        }
        if (!PyBytes_Check(pValue)) {
            throw std::runtime_error("DV2_Vertex: Python function did not return bytes");
        }
        const size_t ret_size = PyBytes_Size(pValue);
        if (ret_size % sizeof(FLOAT) != 0) {
            throw std::runtime_error("DV2_Vertex: returned bytes are not a whole number of floats");
        }
        std::vector<FLOAT> ret(ret_size / sizeof(FLOAT));
        if (ret_size) std::memcpy(ret.data(), PyBytes_AsString(pValue), ret_size);
        return ret;
    }

    std::vector<FLOAT> DV2_Vertex(const std::string& module, const std::string& function,
                                  const std::string& weights,
                                  const std::vector<std::vector<FLOAT> >& input, int n_vertex_rows,
                                  const std::vector<FLOAT>& trad, double q_scale, double q_offset,
                                  const std::string& dtype)
    {
        return dv2_call(module, function, weights, input, n_vertex_rows, trad, q_scale, q_offset, dtype,
                        nullptr, nullptr);
    }

    std::vector<FLOAT> DV2_Vertex(const std::string& module, const std::string& function,
                                  const std::string& weights,
                                  const std::vector<std::vector<FLOAT> >& input, int n_vertex_rows,
                                  const std::vector<FLOAT>& trad, double q_scale, double q_offset,
                                  const std::string& dtype,
                                  const std::vector<std::vector<FLOAT> >& space,
                                  const std::vector<int>& space_cluster_id)
    {
        return dv2_call(module, function, weights, input, n_vertex_rows, trad, q_scale, q_offset, dtype,
                        &space, &space_cluster_id);
    }
}  // namespace WCPPyUtil
