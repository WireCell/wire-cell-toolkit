#ifndef WCPPYUTIL_INTERPRETER_H
#define WCPPYUTIL_INTERPRETER_H

namespace WCPPyUtil {
    /// Start the embedded Python interpreter once per process (thread-safe,
    /// idempotent) and release the GIL afterwards, so that any thread can then
    /// take it with PyGILState_Ensure().  SCN_Vertex, DV2_Load and DV2_Vertex
    /// call this themselves and hold the GIL only for the duration of a call;
    /// other code that calls the Python C API directly must take the GIL the
    /// same way.  If the interpreter was already started by someone else, it is
    /// left as it is (its owner decides who holds the GIL).
    void ensure_interpreter();
}

#endif
