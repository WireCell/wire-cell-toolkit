// Regression test for the libtorch/miniz symbol collision.
//
// WCT vendors miniz (a zip/deflate library) for its .npz / custard I/O.
// libtorch *also* statically bundles its own copy of miniz.  Both used to
// expose the same global C symbols (mz_*, tinfl_*, tdefl_*).  When
// libWireCellUtil.so is loaded into a process before libtorch (as happens in
// wire-cell, which loads Util and then dlopens the WireCellPytorch plugin), the
// dynamic linker bound torch's internal miniz calls to WCT's incompatible copy.
// The two lay out mz_zip_archive differently, so torch::jit::load() dereferenced
// a wild pointer and segfaulted deep inside PyTorchStreamReader::init().
//
// The bug is invisible to test_minimal_torch because that program links only
// libtorch -- there is no WCT miniz in its process to collide.  This test
// deliberately pulls in libWireCellUtil (whose miniz used to be exported) *and*
// libtorch, then loads a model.  Before the fix (WCTMiniz.cmake builds miniz as
// a hidden-visibility static so WCT never exports it) this SIGSEGVs; after, it
// loads cleanly.
//
// Run as an atomic test with no arguments (locates the bundled model), or pass
// a TorchScript model path as argv[1].

#include "WireCellPytorch/Torch.h"

// Reference WireCellUtil so libWireCellUtil.so becomes a load-time dependency of
// this program and its symbols enter the process's global scope -- the very
// condition that used to hijack torch's miniz.
#include "WireCellUtil/Persist.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Locate the bundled TorchScript model.  The runner's working directory has
// varied (build/ vs build/pytorch/), so first derive the model's absolute path
// from this source file's directory (extsmod.ts sits beside it in the source
// tree), then fall back to a few cwd-relative candidates.
static std::string find_extsmod()
{
    std::vector<std::string> candidates;
    {
        const std::string self = __FILE__;         // absolute at compile time
        const auto slash = self.find_last_of('/');
        if (slash != std::string::npos) {
            candidates.push_back(self.substr(0, slash) + "/extsmod.ts");
        }
    }
    candidates.push_back("../pytorch/test/extsmod.ts");     // cwd = build/
    candidates.push_back("../../pytorch/test/extsmod.ts");  // cwd = build/pytorch/
    candidates.push_back("pytorch/test/extsmod.ts");        // cwd = source top
    for (const auto& path : candidates) {
        if (std::ifstream(path).good()) {
            return path;
        }
    }
    return "";
}

int main(int argc, char* argv[])
{
    std::string ts = argc > 1 ? argv[1] : find_extsmod();
    if (ts.empty()) {
        // No model to load: nothing to exercise.  Skip rather than fail so the
        // suite stays green in checkouts without the bundled model.
        std::cerr << "test_torch_miniz_collision: model not found, skipping\n";
        return 0;
    }

    // Force libWireCellUtil to actually be needed at link/load time by touching
    // one of its symbols.  (The result is irrelevant.)
    const bool exists = WireCell::Persist::exists(ts);
    std::cout << "using model " << ts << " (exists=" << exists << ")\n";

    // The load below is what used to segfault.  A genuine c10::Error (bad/missing
    // model) is a *clean* failure, not the crash this test guards against.
    try {
        torch::jit::script::Module mod = torch::jit::load(ts);
        std::cout << "loaded model without miniz collision\n";
    }
    catch (const c10::Error& e) {
        std::cerr << "load raised (not a crash): " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
