# WCT Dev Area Configuration with wcwc

## Overview

This documents how to configure a Wire-Cell Toolkit (WCT) development area using `wcwc`
(the Wire-Cell spack environment manager), including how to set up a view with TMVA-enabled ROOT.

---

## Standard configure command

```bash
ROOTSYS=$(root-config --prefix) ./wcb configure --prefix=$PREFIX \
     --with-cuda=no --with-libtorch=no --with-root=yes \
     --boost-mt --boost-libs=$PREFIX/lib --boost-include=$PREFIX/include \
     --with-jsonnet-libs=gojsonnet --build-debug="-O2 -ggdb3"
```

**Why `--with-root=yes` instead of a path:** the `wcb` root detection gate (`with_p('root')`)
only activates when some `--with-root*` flag is present. `--with-root=yes` signals "detect ROOT
automatically" without hardcoding a path. The actual location is resolved from `ROOTSYS`.

**Why `ROOTSYS=$(root-config --prefix)`:** `wcb` searches `ROOTSYS/bin` for `root-config`.
Setting it inline from the view's `root-config` ensures the build uses exactly the ROOT that
the view provides.

---

## Updating the view to a different spack spec (e.g. to get ROOT+TMVA)

The `local/` directory is a spack view — a flat tree of symlinks into spack-managed packages.
When you need a different set of packages (e.g. ROOT built with TMVA support), recreate the view.

### 1. Find the spack hash to use

```bash
# List all ROOT builds with TMVA support
wcwc spack find -lv -p root+tmva

# Find a WCT build and inspect its full dependency tree
wcwc spack find -lv -p wire-cell-toolkit
wcwc spack spec /<hash>        # confirm root has +tmva+tmva-cpu in the spec
```

### 2. Back up the current view

```bash
cd <dev-area>        # e.g. ~/work/scratch_wcgpu1/toolkit-dev
mv local local-old
```

### 3. Recreate the view from the target spack spec

```bash
wcwc view -S wirecell -s "/<hash>" local
```

`-S wirecell` selects the wirecell view profile; `-s "/<hash>"` pins the root spec.

### 4. Remove the bundled WCT from the view

The view includes the spack-installed WireCell libraries. Remove them so your own build
(under `build/`) is the only WireCell on the path:

```bash
rm -r local/lib/libWireCell* local/include/WireCell*
```

### 5. Verify ROOT has TMVA

```bash
local/bin/root-config --features | tr ' ' '\n' | grep tmva
# Should print: tmva  tmva-cpu
```

### 6. Reconfigure and rebuild

```bash
ROOTSYS=$(root-config --prefix) ./wcb configure --prefix=$PREFIX \
     --with-cuda=no --with-libtorch=no --with-root=yes \
     --boost-mt --boost-libs=$PREFIX/lib --boost-include=$PREFIX/include \
     --with-jsonnet-libs=gojsonnet --build-debug="-O2 -ggdb3"

wcbuild
```

### Rollback

If something goes wrong, restore the old view:

```bash
rm -rf local && mv local-old local
```

---

## Notes

- The view symlinks all runtime dependencies of the chosen spack spec. After `wcwc view`,
  `local/bin/root-config` resolves to whichever ROOT that spec was built against.
- The `wcwc spack spec /<hash>` command prints the full dependency tree including the ROOT
  variant flags (`+tmva`, `+tmva-cpu`, `+opengl`, etc.). Use it to confirm before switching.
- Two ROOT builds with TMVA are available (as of 2026-04):
  - `5umcjrn`: `root@6.32.02 +tmva+tmva-cpu ~opengl` — used by WCT 0.34.2 spack build `/2u5bdsx`
  - `kl5ziao`: `root@6.32.02 +tmva+tmva-cpu +opengl`
