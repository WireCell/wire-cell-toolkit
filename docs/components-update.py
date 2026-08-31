#!/usr/bin/env python3
"""
components-update.py — maintenance tool for docs/components.md and docs/components/*.md

Run from the toolkit root directory.

Subcommands
-----------
  inventory   Show new, removed, and source-modified components
  regen       Rebuild docs/components.md from existing docs/components/*.md files
  write FILE  Parse a component block from FILE (or stdin if FILE is -) and write
              docs/components/NAME.md, then update the docs/components.md row

The 'write' subcommand consumes the structured block format described in
docs/components-update.md and produced by the LLM update prompt.
"""

import argparse
import glob
import os
import re
import subprocess
import sys

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
TOOLKIT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(TOOLKIT, "docs")
COMP_DIR = os.path.join(DOCS, "components")
COMP_MD = os.path.join(DOCS, "components.md")

# ---------------------------------------------------------------------------
# Interface → (category, input, output) table
# Derived from iface/inc/WireCellIface/*.h
# ---------------------------------------------------------------------------
IFACE_MAP = {
    "IFrameFilter":         ("function",   "IFrame",           "IFrame"),
    "IFrameSource":         ("source",     "(none)",           "IFrame"),
    "IFrameSink":           ("sink",       "IFrame",           "(none)"),
    "IFrameFanin":          ("fanin",      "vector<IFrame>",   "IFrame"),
    "IFrameFanout":         ("fanout",     "IFrame",           "vector<IFrame>"),
    "IFrameJoiner":         ("join",       "(IFrame, IFrame)", "IFrame"),
    "IFrameMerge":          ("hydra",      "vector<IFrame>",   "IFrame"),
    "IFrameSlicer":         ("queuedout",  "IFrame",           "ISlice"),
    "IFrameSlices":         ("function",   "IFrame",           "ISliceSet"),
    "IFrameSplitter":       ("split",      "IFrame",           "(IFrame, IFrame)"),
    "IFrameTensorSet":      ("function",   "IFrame",           "ITensorSet"),
    "IDepoFilter":          ("function",   "IDepo",            "IDepo"),
    "IDepoSource":          ("source",     "(none)",           "IDepo"),
    "IDepoSink":            ("sink",       "IDepo",            "(none)"),
    "IDepoFanout":          ("fanout",     "IDepo",            "vector<IDepo>"),
    "IDepoMerger":          ("hydra",      "(IDepo, IDepo)",   "IDepo"),
    "IDepoCollector":       ("queuedout",  "IDepo",            "IDepoSet"),
    "IDepoFramer":          ("function",   "IDepoSet",         "IFrame"),
    "IDepoSetFilter":       ("function",   "IDepoSet",         "IDepoSet"),
    "IDepoSetSource":       ("source",     "(none)",           "IDepoSet"),
    "IDepoSetSink":         ("sink",       "IDepoSet",         "(none)"),
    "IDepoSetFanin":        ("fanin",      "vector<IDepoSet>", "IDepoSet"),
    "IDepoSetFanout":       ("fanout",     "IDepoSet",         "vector<IDepoSet>"),
    "IDepos2DeposOrFrame":  ("hydra",      "IDepoSet",         "IDepoSet or IFrame"),
    "IClusterFilter":       ("function",   "ICluster",         "ICluster"),
    "IClusterSource":       ("source",     "(none)",           "ICluster"),
    "IClusterSink":         ("sink",       "ICluster",         "(none)"),
    "IClusterFanin":        ("fanin",      "vector<ICluster>", "ICluster"),
    "IClusterFanout":       ("fanout",     "ICluster",         "vector<ICluster>"),
    "IClusterFramer":       ("function",   "ICluster",         "IFrame"),
    "IClusterTensorSet":    ("function",   "ICluster",         "ITensorSet"),
    "IClusterFaninTensorSet": ("fanin",    "vector<ICluster>", "ITensorSet"),
    "IBlobSetProcessor":    ("function",   "IBlobSet",         "IBlobSet"),
    "IBlobSetSink":         ("sink",       "IBlobSet",         "(none)"),
    "IBlobSetSource":       ("source",     "(none)",           "IBlobSet"),
    "IBlobSetFanin":        ("fanin",      "vector<IBlobSet>", "IBlobSet"),
    "IBlobSetFanout":       ("fanout",     "IBlobSet",         "vector<IBlobSet>"),
    "IBlobSetFramer":       ("function",   "IBlobSet",         "IFrame"),
    "IBlobDeclustering":    ("function",   "ICluster",         "IBlobSet"),
    "IBlobDepoFill":        ("join",       "(ICluster, IDepoSet)", "ICluster"),
    "IBlobSampling":        ("function",   "IBlobSet",         "ITensorSet"),
    "IBlobTensoring":       ("queuedout",  "IBlobSet",         "ITensorSet"),
    "ITensorSetFilter":     ("function",   "ITensorSet",       "ITensorSet"),
    "ITensorSetSource":     ("source",     "(none)",           "ITensorSet"),
    "ITensorSetSink":       ("sink",       "ITensorSet",       "(none)"),
    "ITensorSetFanin":      ("fanin",      "vector<ITensorSet>", "ITensorSet"),
    "ITensorSetFrame":      ("function",   "ITensorSet",       "IFrame"),
    "ITensorSetCluster":    ("function",   "ITensorSet",       "ICluster"),
    "ITensorSetUnpacker":   ("fanout",     "ITensorSet",       "vector<ITensor>"),
    "ITensorPacker":        ("fanin",      "vector<ITensor>",  "ITensorSet"),
    "ISliceFanout":         ("fanout",     "ISlice",           "vector<ISlice>"),
    "ISliceFrameSink":      ("sink",       "(ISlice, IFrame)", "(none)"),
    "ISliceStriper":        ("function",   "ISlice",           "IStripeSet"),
    "ITiling":              ("function",   "ISlice",           "IBlobSet"),
    "IWireGenerator":       ("function",   "IWireParameters",  "IWire::vector"),
    "IWireSource":          ("source",     "(none)",           "IWire::vector"),
    "IWireSummarizer":      ("function",   "IWire::vector",    "IWireSummary"),
    "IClustering":          ("queuedout",  "IBlobSet",         "ICluster"),
    "IDuctor":              ("queuedout",  "IDepo",            "IFrame"),
    "IDiffuser":            ("queuedout",  "IDepo",            "IDiffusion"),
    "IDrifter":             ("queuedout",  "IDepo",            "IDepo"),
    "IPointFieldSink":      ("sink",       "IPointField",      "(none)"),
    "IScalarFieldSink":     ("sink",       "IScalarField",     "(none)"),
}

# Interfaces that indicate INode participation
INODE_IFACES = set(IFACE_MAP.keys())

# ---------------------------------------------------------------------------
# Source inventory
# ---------------------------------------------------------------------------

def find_factory_components():
    """Return dict name -> {name, concrete, ifaces, srcfile} for all
    WIRECELL_FACTORY registrations that include an INode-derived interface."""
    pattern = re.compile(r'WIRECELL_FACTORY\s*\((.+?)\)', re.DOTALL)
    result = {}
    for srcfile in sorted(glob.glob(os.path.join(TOOLKIT, '**/*.cxx'), recursive=True)):
        try:
            content = open(srcfile).read()
        except Exception:
            continue
        if 'WIRECELL_FACTORY' not in content:
            continue
        for m in pattern.findall(content):
            args = [a.strip() for a in m.split(',')]
            if len(args) < 2:
                continue
            name = args[0]
            concrete = args[1]
            ifaces = args[2:]
            # Check for INode-derived interface (strip WireCell:: namespace)
            bare_ifaces = [i.split('::')[-1] for i in ifaces]
            inode_iface = next((i for i in bare_ifaces if i in INODE_IFACES), None)
            if inode_iface is None:
                continue
            result[name] = {
                'name': name,
                'concrete': concrete,
                'ifaces': ifaces,
                'inode_iface': inode_iface,
                'srcfile': srcfile,
            }
    return result


def infer_properties(inode_iface):
    """Return (category, input, output) from the interface name."""
    return IFACE_MAP.get(inode_iface, ('unknown', 'unknown', 'unknown'))

# ---------------------------------------------------------------------------
# docs/components/*.md parsing
# ---------------------------------------------------------------------------

def parse_component_page(path):
    """Parse a NAME.md file.  Returns dict with keys:
    name, description, concrete, category, iface, input, output,
    configurable, params (list of (param, desc) tuples)."""
    text = open(path).read()
    lines = text.split('\n')

    comp = {'params': []}

    # Title = first # heading
    for line in lines:
        m = re.match(r'^#\s+(.+)', line)
        if m:
            comp['name'] = m.group(1).strip()
            break

    # Description = first non-empty paragraph after the title
    after_title = False
    desc_lines = []
    for line in lines:
        if re.match(r'^#\s+', line):
            after_title = True
            continue
        if not after_title:
            continue
        if re.match(r'^##', line):
            break
        if line.strip():
            desc_lines.append(line.strip())
        elif desc_lines:
            break
    comp['description'] = ' '.join(desc_lines)

    # Node Properties table
    prop_map = {
        'Factory name': 'name',
        'Concrete class': 'concrete',
        'Node category': 'category',
        'Primary interface': 'iface',
        'Input type(s)': 'input',
        'Output type(s)': 'output',
        'Configurable': 'configurable',
    }
    for line in lines:
        for label, key in prop_map.items():
            if f'| {label} |' in line:
                val = line.split('|')[2].strip().strip('`')
                comp[key] = val

    # Configuration Parameters table
    in_params = False
    for line in lines:
        if '## Configuration Parameters' in line:
            in_params = True
            continue
        if in_params and line.startswith('## '):
            break
        if in_params and line.startswith('| `'):
            parts = line.split('|')
            if len(parts) >= 3:
                param = parts[1].strip().strip('`')
                desc = parts[2].strip()
                comp['params'].append((param, desc))

    return comp


def read_all_component_pages():
    """Return dict name -> parsed component dict for all existing NAME.md files."""
    result = {}
    for path in sorted(glob.glob(os.path.join(COMP_DIR, '*.md'))):
        try:
            comp = parse_component_page(path)
            if comp.get('name'):
                result[comp['name']] = comp
        except Exception as e:
            print(f"Warning: could not parse {path}: {e}", file=sys.stderr)
    return result

# ---------------------------------------------------------------------------
# Writing docs
# ---------------------------------------------------------------------------

def write_component_page(comp):
    """Write docs/components/NAME.md from a component dict."""
    name = comp['name']
    configurable = str(comp.get('configurable', 'no')).lower() == 'yes'
    params = comp.get('params', [])
    real_params = [(p, d) for p, d in params
                   if p and p.lower() not in ('none', '(none)')]

    lines = [
        f"# {name}\n",
        f"{comp.get('description', '')}\n",
        "## Node Properties\n",
        "| Property | Value |",
        "| --- | --- |",
        f"| Factory name | `{name}` |",
        f"| Concrete class | `{comp.get('concrete', '')}` |",
        f"| Node category | {comp.get('category', '')} |",
        f"| Primary interface | `{comp.get('iface', '')}` |",
        f"| Input type(s) | `{comp.get('input', '')}` |",
        f"| Output type(s) | `{comp.get('output', '')}` |",
        f"| Configurable | {'yes' if configurable else 'no'} |",
        "",
    ]

    if configurable and real_params:
        lines += [
            "## Configuration Parameters\n",
            "| Parameter | Description |",
            "| --- | --- |",
        ]
        for param, desc in real_params:
            lines.append(f"| `{param}` | {desc} |")
        lines.append("")
    elif configurable:
        lines += [
            "## Configuration Parameters\n",
            "This component is configurable but has no documented parameters.\n",
        ]

    os.makedirs(COMP_DIR, exist_ok=True)
    outpath = os.path.join(COMP_DIR, f"{name}.md")
    with open(outpath, 'w') as f:
        f.write('\n'.join(lines))
    return outpath


def regen_summary_table(components):
    """Rewrite docs/components.md from a dict of parsed component pages."""

    def pkg_key(comp):
        c = comp.get('concrete', '')
        order = ['Gen', 'Sig::', 'SigProc', 'Img', 'Clus', 'Aux', 'Sio', 'Root', 'Hio', 'Pytorch', 'Zio']
        for i, token in enumerate(order):
            if token in c:
                return (i, comp['name'])
        return (99, comp['name'])

    sorted_comps = sorted(components.values(), key=pkg_key)

    lines = [
        "# Wire-Cell Toolkit Data-Flow Components\n",
        "This file lists all components that participate as nodes in the WCT data-flow graph.",
        "Each component is registered with `WIRECELL_FACTORY` and implements an `INode`-derived interface.\n",
        "| Name | Concrete Class | Category | Interface | Input | Output | Description |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]

    for comp in sorted_comps:
        name = comp['name']
        concrete = comp.get('concrete', '')
        cat = comp.get('category', '')
        iface = comp.get('iface', '')
        inp = comp.get('input', '').replace('|', '\\|')
        out = comp.get('output', '').replace('|', '\\|')
        desc = comp.get('description', '')
        if len(desc) > 100:
            desc = desc[:97] + '...'
        lines.append(
            f"| [{name}](components/{name}.md) | `{concrete}` | {cat}"
            f" | `{iface}` | `{inp}` | `{out}` | {desc} |"
        )

    lines.append("")
    with open(COMP_MD, 'w') as f:
        f.write('\n'.join(lines))
    print(f"Wrote {COMP_MD} ({len(sorted_comps)} entries)")

# ---------------------------------------------------------------------------
# Block format parser  (input for 'write' subcommand)
# ---------------------------------------------------------------------------

def parse_block(text):
    """Parse the LLM output block format into a component dict.

    Expected format (fields in any order, terminated by --- or EOF):

        NAME: Foo
        CONCRETE: WireCell::Bar::Foo
        CATEGORY: function
        INPUT: IFrame
        OUTPUT: IFrame
        IFACE: IFrameFilter
        DESCRIPTION: One line description.
        CONFIGURABLE: yes
        CONFIG_PARAMS:
        param1: description of param1
        param2: description of param2
    """
    comp = {'params': []}
    in_config = False

    for line in text.split('\n'):
        line = line.rstrip()
        if line == '---':
            break
        elif line.startswith('NAME:'):
            comp['name'] = line[5:].strip()
        elif line.startswith('CONCRETE:'):
            comp['concrete'] = line[9:].strip()
        elif line.startswith('CATEGORY:'):
            comp['category'] = line[9:].strip()
        elif line.startswith('INPUT:'):
            comp['input'] = line[6:].strip()
        elif line.startswith('OUTPUT:'):
            comp['output'] = line[7:].strip()
        elif line.startswith('IFACE:'):
            comp['iface'] = line[6:].strip()
        elif line.startswith('DESCRIPTION:'):
            comp['description'] = line[12:].strip()
        elif line.startswith('CONFIGURABLE:'):
            comp['configurable'] = line[13:].strip()
        elif line.startswith('CONFIG_PARAMS:'):
            in_config = True
            rest = line[14:].strip()
            if rest and rest.lower() not in ('none', '(none)'):
                if ':' in rest:
                    p, d = rest.split(':', 1)
                    comp['params'].append((p.strip(), d.strip()))
        elif in_config and line:
            if ':' in line:
                p, d = line.split(':', 1)
                comp['params'].append((p.strip(), d.strip()))

    return comp

# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_inventory(args):
    """Compare source against docs and report differences."""
    source = find_factory_components()
    documented = read_all_component_pages()

    src_names = set(source.keys())
    doc_names = set(documented.keys())

    new = sorted(src_names - doc_names)
    removed = sorted(doc_names - src_names)
    common = src_names & doc_names

    if new:
        print(f"\nNEW ({len(new)}) — need docs/components/NAME.md and a row in components.md:")
        for n in new:
            info = source[n]
            print(f"  {n}")
            print(f"    src:   {info['srcfile']}")
            print(f"    iface: {info['inode_iface']}")

    if removed:
        print(f"\nREMOVED ({len(removed)}) — delete docs/components/NAME.md and its row:")
        for n in removed:
            print(f"  {n}")

    # Check for interface changes in common components
    changed = []
    for n in sorted(common):
        src = source[n]
        doc = documented[n]
        src_iface = src['inode_iface']
        doc_iface = doc.get('iface', '')
        src_concrete = src['concrete'].split('::')[-1]  # just class name
        doc_concrete = doc.get('concrete', '').split('::')[-1]
        if src_iface != doc_iface or src_concrete != doc_concrete:
            changed.append((n, src, doc))

    if changed:
        print(f"\nSTRUCTURALLY CHANGED ({len(changed)}) — interface or class name differs:")
        for n, src, doc in changed:
            print(f"  {n}")
            if src['inode_iface'] != doc.get('iface'):
                print(f"    iface: {doc.get('iface')} -> {src['inode_iface']}")
            if src['concrete'] != doc.get('concrete'):
                print(f"    concrete: {doc.get('concrete')} -> {src['concrete']}")

    # Show recently git-modified source files for common components
    if not args.no_git:
        try:
            result = subprocess.run(
                ['git', 'diff', '--name-only', 'HEAD~1'],
                capture_output=True, text=True, cwd=TOOLKIT
            )
            changed_files = set(result.stdout.split())
            src_modified = []
            for n in sorted(common):
                rel = os.path.relpath(source[n]['srcfile'], TOOLKIT)
                if rel in changed_files:
                    src_modified.append((n, source[n]['srcfile']))
            if src_modified:
                print(f"\nSOURCE MODIFIED SINCE HEAD~1 ({len(src_modified)}) — review and update docs:")
                for n, f in src_modified:
                    print(f"  {n}  ({f})")
        except Exception:
            pass  # git not available or no prior commit

    if not new and not removed and not changed:
        print("No structural changes detected.")


def cmd_regen(args):
    """Rebuild docs/components.md from all existing NAME.md files."""
    components = read_all_component_pages()
    print(f"Read {len(components)} component pages.")
    regen_summary_table(components)


def cmd_write(args):
    """Parse a block file and write the NAME.md and update components.md."""
    if args.file == '-':
        text = sys.stdin.read()
    else:
        text = open(args.file).read()

    comp = parse_block(text)
    if not comp.get('name'):
        print("Error: no NAME field found in block.", file=sys.stderr)
        sys.exit(1)

    # Fill in category/input/output from IFACE_MAP if not specified
    iface = comp.get('iface', '')
    if iface and ('category' not in comp or 'input' not in comp or 'output' not in comp):
        cat, inp, out = infer_properties(iface)
        comp.setdefault('category', cat)
        comp.setdefault('input', inp)
        comp.setdefault('output', out)

    outpath = write_component_page(comp)
    print(f"Wrote {outpath}")

    # Update the row in components.md
    components = read_all_component_pages()
    regen_summary_table(components)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest='cmd', required=True)

    p_inv = sub.add_parser('inventory', help='Show new/removed/changed components')
    p_inv.add_argument('--no-git', action='store_true',
                       help='Skip git-based source-change detection')

    sub.add_parser('regen', help='Rebuild docs/components.md from existing NAME.md files')

    p_write = sub.add_parser('write',
                             help='Parse a block file and write NAME.md + update table')
    p_write.add_argument('file', help='Block file path, or - to read from stdin')

    args = parser.parse_args()
    {'inventory': cmd_inventory, 'regen': cmd_regen, 'write': cmd_write}[args.cmd](args)


if __name__ == '__main__':
    main()
