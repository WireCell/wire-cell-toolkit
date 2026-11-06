#include <Python.h>

#include "WCPPyUtil/Interpreter.h"

#include <mutex>

void WCPPyUtil::ensure_interpreter()
{
    static std::once_flag once;
    std::call_once(once, [] {
        if (Py_IsInitialized()) return;   // started elsewhere: leave its GIL as it is
        Py_Initialize();                  // never finalised (torch / sparseconvnet crash in Py_Finalize)
        PyEval_SaveThread();              // release the GIL the initialising thread now holds
    });
}
