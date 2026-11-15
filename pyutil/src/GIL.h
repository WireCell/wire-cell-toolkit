// Private to WCPPyUtil: the embedded interpreter's GIL for the scope of one
// call into Python (from any thread).  Declare it before any PyObject handling
// (e.g. a PyRefs holder) so that references are released while it is held.
#ifndef WCPPYUTIL_GIL_H
#define WCPPYUTIL_GIL_H

#include <Python.h>

#include "WCPPyUtil/Interpreter.h"

namespace WCPPyUtil {
    class GIL {
      public:
        GIL()
        {
            ensure_interpreter();
            m_state = PyGILState_Ensure();
        }
        ~GIL() { PyGILState_Release(m_state); }
        GIL(const GIL&) = delete;
        GIL& operator=(const GIL&) = delete;

      private:
        PyGILState_STATE m_state;
    };
}

#endif
