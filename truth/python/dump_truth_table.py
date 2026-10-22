#!/usr/bin/env python3
"""dump_truth_table.py -- Stage A truth dumper: the MC particle tables of one art file -> truth_<ident>.npz per event.

Generalised from wcp-porting-img fdvd_sim/stageA/dump_truth.py (gallery, FD-VD solar), for every detector
whose Stage A runs in LArSoft (FD-HD, PD-HD, PD-VD, FD-VD with dunesw; SBND reco1 with sbndcode): the
MCTruth products of every generator label, every MCParticle with Process(), momentum and the generator
join through the MCTruth<->MCParticle Assns (descendants inherit their primary's label through the mother
chain), the SimEnergyDeposit ids (TrackID, OrigTrackID per deposit, in vector order = the depo id of
wclsSimDepoSetSource index mode, owner decision O3) and the visible electrons per trackid.

usage (inside the experiment's LArSoft environment, gallery on the ROOT path):
  python3 -I dump_truth_table.py <art.root> <outdir> [--entries 0,2-5] [--ident entry|event]
                                 [--sed-tag IonAndScint[:instance]] [--mcp-tag largeant]
                                 [--detector NAME] [--id-mode index|track] [--sce none|...] [--meta k=v ...]
Output: <outdir>/truth_<ident>.npz (ident = entry index, or the art event number with --ident event) and
<outdir>/truth_<ident>.txt (one summary line); existing outputs are not overwritten (add-only).

Columns and units: truth_table.py (SBND attacher conventions, cm / ns / GeV).  nu_row of a particle is the
truth_nu row of its generator interaction; a non-neutrino generator (cosmics, single particles,
radiologicals) gets one truth_nu row per label with pdg 0 and nu_idx -1.  parent_trackid equals
mother_trackid at dump time (every particle is kept; PFMerger recomputes the kept parent).  pf_node, cls,
merged_into are -1 placeholders filled by PFMerger; nvis_q is the IonAndScint electron sum of the trackid.
"""
import argparse, collections, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir
import truth_table as T  # noqa: E402


def parse_entries(spec, n):
    if not spec:
        return list(range(n))
    out = []
    for part in spec.split(","):
        if "-" in part:
            a, b = part.split("-"); out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return [i for i in out if 0 <= i < n]


def branch_labels(fn, prefix):
    """(label, instance, process) of every Events branch named <prefix><label>_<instance>_<process>."""
    import ROOT
    f = ROOT.TFile.Open(fn); t = f.Get("Events")
    out = []
    for br in t.GetListOfBranches():
        n = br.GetName()
        if n.startswith(prefix):
            parts = n[len(prefix):].rstrip(".").split("_")
            if len(parts) >= 3:
                out.append((parts[0], parts[1], parts[2]))
    f.Close()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("artfile"); ap.add_argument("outdir")
    ap.add_argument("--entries", default="")
    ap.add_argument("--ident", default="entry", choices=["entry", "event"])
    ap.add_argument("--sed-tag", default="IonAndScint", help="label[:instance] of the SimEnergyDeposits the depos came from")
    ap.add_argument("--mcp-tag", default="", help="MCParticle label (default: the first simb::MCParticles branch)")
    ap.add_argument("--detector", default="")
    ap.add_argument("--id-mode", default="index", choices=["index", "track"], help="wclsSimDepoSetSource id_is_track false|true (recorded)")
    ap.add_argument("--sce", default="none", help="SCE status of the Stage A deposits (recorded; owner decision O4)")
    ap.add_argument("--meta", nargs="*", default=[], help="extra k=v provenance")
    ap.add_argument("--prefix", default="truth")
    a = ap.parse_args()

    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gSystem.Load("libgallery")
    ROOT.gInterpreter.ProcessLine('#include "gallery/Event.h"')
    ROOT.gInterpreter.ProcessLine('#include "canvas/Persistency/Common/Assns.h"')
    ROOT.gInterpreter.ProcessLine('#include "nusimdata/SimulationBase/MCTruth.h"')
    ROOT.gInterpreter.ProcessLine('#include "nusimdata/SimulationBase/MCParticle.h"')
    ROOT.gInterpreter.ProcessLine('#include "lardataobj/Simulation/SimEnergyDeposit.h"')

    os.makedirs(a.outdir, exist_ok=True)
    mct_labels = branch_labels(a.artfile, "simb::MCTruths_")
    mcp_labels = branch_labels(a.artfile, "simb::MCParticles_")
    sed_labels = branch_labels(a.artfile, "sim::SimEnergyDeposits_")
    if not mcp_labels:
        sys.exit("no simb::MCParticles branch in " + a.artfile)
    mcp_lab = a.mcp_tag or mcp_labels[0][0]
    sed_lab, sed_inst = (a.sed_tag.split(":", 1) + [""])[:2]
    have_sed = any(l == sed_lab and i == sed_inst for l, i, p in sed_labels)
    if not have_sed:
        print("WARNING: no SimEnergyDeposits %s:%s in the file (have %s); sed tables will be empty"
              % (sed_lab, sed_inst, sed_labels))
    have_anc = True

    ev = ROOT.gallery.Event(ROOT.std.vector("string")(1, a.artfile))
    nentries = ev.numberOfEventsInFile()
    wanted = set(parse_entries(a.entries, nentries))
    meta_extra = dict(kv.split("=", 1) for kv in a.meta)

    def handle(cls, label, instance=""):
        return ev.getValidHandle[cls](ROOT.art.InputTag(label, instance))

    entry = -1
    while not ev.atEnd():
        entry += 1
        if entry not in wanted:
            ev.next(); continue
        aux = ev.eventAuxiliary()
        run, sub, evt = aux.run(), aux.subRun(), aux.event()
        ident = entry if a.ident == "entry" else evt
        out = os.path.join(a.outdir, "%s_%d.npz" % (a.prefix, ident))
        if os.path.exists(out):
            print("entry %d (r%d s%d e%d): %s exists, skipped" % (entry, run, sub, evt, out)); ev.next(); continue

        # ---- truth_nu: one row per MCTruth with a neutrino; one per non-neutrino generator label ----
        gen_labels = [l for l, i, p in mct_labels]
        pid_to_row = {}        # (product id, key) -> truth_nu row
        label_row = {}         # label idx -> truth_nu row of a non-neutrino label
        nu_rows = []
        nu_idx = 0
        for li, (lab, inst, proc) in enumerate(mct_labels):
            h = handle("vector<simb::MCTruth>", lab, inst)
            pid = h.id().value()
            for k in range(h.product().size()):
                mct = h.product()[k]
                row = np.zeros(len(T.NU_COLUMNS))
                if mct.NeutrinoSet():
                    nu = mct.GetNeutrino(); p = nu.Nu()
                    row[T.nu_col("nu_idx")] = nu_idx; nu_idx += 1
                    row[T.nu_col("pdg")] = p.PdgCode(); row[T.nu_col("ccnc")] = nu.CCNC(); row[T.nu_col("mode")] = nu.Mode()
                    row[T.nu_col("int_type")] = nu.InteractionType(); row[T.nu_col("flavor")] = abs(p.PdgCode())
                    row[T.nu_col("E")] = p.E(); row[T.nu_col("vtx_x")] = p.Vx(); row[T.nu_col("vtx_y")] = p.Vy()
                    row[T.nu_col("vtx_z")] = p.Vz(); row[T.nu_col("t")] = p.T()
                    pid_to_row[(pid, k)] = len(nu_rows); nu_rows.append(row)
                else:
                    if li not in label_row:
                        row[T.nu_col("nu_idx")] = -1; row[T.nu_col("pdg")] = 0
                        if mct.NParticles() > 0:
                            p = mct.GetParticle(0)
                            row[T.nu_col("vtx_x")] = p.Vx(); row[T.nu_col("vtx_y")] = p.Vy(); row[T.nu_col("vtx_z")] = p.Vz()
                            row[T.nu_col("t")] = p.T()
                        label_row[li] = len(nu_rows); nu_rows.append(row)
                    pid_to_row[(pid, k)] = label_row[li]
        row_label = {}
        for li, (lab, inst, proc) in enumerate(mct_labels):
            h = handle("vector<simb::MCTruth>", lab, inst); pid = h.id().value()
            for k in range(h.product().size()):
                row_label[pid_to_row[(pid, k)]] = li

        # ---- truth_pf: every MCParticle ----
        hp = handle("vector<simb::MCParticle>", mcp_lab)
        parts = hp.product(); n = parts.size()
        pf = np.zeros((n, len(T.PF_COLUMNS)))
        row_of = {}
        for i in range(n):
            p = parts[i]
            nt = p.NumberTrajectoryPoints()
            tid = p.TrackId(); row_of[tid] = i
            pf[i, T.pf_col("trackid")] = tid; pf[i, T.pf_col("mother_trackid")] = p.Mother(); pf[i, T.pf_col("parent_trackid")] = p.Mother()
            pf[i, T.pf_col("pdg")] = p.PdgCode(); pf[i, T.pf_col("process")] = T.process_code(str(p.Process()))
            pf[i, T.pf_col("E")] = p.E(0); pf[i, T.pf_col("KE")] = p.E(0) - p.Mass()
            pf[i, T.pf_col("start_x")] = p.Vx(0); pf[i, T.pf_col("start_y")] = p.Vy(0); pf[i, T.pf_col("start_z")] = p.Vz(0); pf[i, T.pf_col("start_t")] = p.T(0)
            pf[i, T.pf_col("start_px")] = p.Px(0); pf[i, T.pf_col("start_py")] = p.Py(0); pf[i, T.pf_col("start_pz")] = p.Pz(0)
            j = nt - 1 if nt > 0 else 0
            pf[i, T.pf_col("end_x")] = p.Vx(j); pf[i, T.pf_col("end_y")] = p.Vy(j); pf[i, T.pf_col("end_z")] = p.Vz(j); pf[i, T.pf_col("end_t")] = p.T(j)
            pf[i, T.pf_col("orig_trackid")] = tid
            pf[i, T.pf_col("nu_row")] = -1; pf[i, T.pf_col("gen_label")] = -1
            pf[i, T.pf_col("pf_node")] = -1; pf[i, T.pf_col("cls")] = T.CLS["unknown"]; pf[i, T.pf_col("merged_into")] = -1
        # generator join: Assns<MCTruth, MCParticle, GeneratedParticleInfo> then the mother chain
        tid_row = {}
        try:
            A = handle("art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo>", mcp_lab).product()
            for j in range(A.size()):
                pr = A.at(j)
                r = pid_to_row.get((pr.first.id().value(), pr.first.key()), -1)
                if pr.second.key() < n:
                    tid_row[parts[pr.second.key()].TrackId()] = r
        except Exception as e:
            print("no MCTruth<->MCParticle Assns (%s): generator join by the mother chain only" % e)
        sys.setrecursionlimit(100000)
        mother = {int(pf[i, T.pf_col("trackid")]): int(pf[i, T.pf_col("mother_trackid")]) for i in range(n)}

        def nu_row_of(tid, depth=0):
            if tid in tid_row: return tid_row[tid]
            m = mother.get(tid, 0)
            r = nu_row_of(m, depth + 1) if (m != 0 and depth < 100000) else -1
            tid_row[tid] = r
            return r
        for i in range(n):
            r = nu_row_of(int(pf[i, T.pf_col("trackid")]))
            pf[i, T.pf_col("nu_row")] = r
            pf[i, T.pf_col("gen_label")] = row_label.get(r, -1)

        # ancestry map (dropped descendants -> kept ancestor), when the G4 stage wrote one
        anc = []
        if have_anc:
            try:
                am = handle("sim::ParticleAncestryMap", mcp_lab).product().GetMap()
                for kv in am:
                    for d in kv.second:
                        anc.append((int(d), int(kv.first)))
            except Exception as e:
                have_anc = False
                print("no ancestry map:", e)
        anc_dict = dict(anc)

        # ---- SimEnergyDeposits: ids per deposit and electrons per trackid ----
        sed_ids = np.zeros((0, 2), np.int32); sed_trackq = np.zeros((0, 2)); sed_energy = np.zeros(0, np.float32)
        edep_row = collections.defaultdict(float)
        if have_sed:
            seds = handle("vector<sim::SimEnergyDeposit>", sed_lab, sed_inst).product()
            ns = seds.size()
            sed_ids = np.zeros((ns, 2), np.int32); sed_energy = np.zeros(ns, np.float32)
            tq = collections.defaultdict(float)
            for k in range(ns):
                s = seds[k]
                tid = int(s.TrackID())
                try:
                    otid = int(s.OrigTrackID())
                except Exception:
                    otid = tid
                sed_ids[k, 0] = tid; sed_ids[k, 1] = otid; sed_energy[k] = s.Energy()
                tq[tid] += s.NumElectrons()
                kept = abs(tid)
                if kept not in row_of and kept in anc_dict: kept = anc_dict[kept]
                r = tid_row.get(kept, -1)
                edep_row[r] += s.Energy()
            sed_trackq = np.array(sorted(tq.items()), dtype=np.float64).reshape(-1, 2)
            for i in range(n):
                pf[i, T.pf_col("nvis_q")] = tq.get(int(pf[i, T.pf_col("trackid")]), 0.0)
        nu = np.array(nu_rows, dtype=np.float64).reshape(-1, len(T.NU_COLUMNS))
        for r, e in edep_row.items():
            if 0 <= r < len(nu): nu[r, T.nu_col("edep")] = e * 1e-3   # MeV -> GeV
        meta = {"detector": a.detector, "artfile": os.path.abspath(a.artfile), "entry": entry,
                "run": run, "subrun": sub, "event": evt, "mcp_tag": mcp_lab, "sed_tag": a.sed_tag,
                "id_mode": a.id_mode, "sce": a.sce, "have_sed": bool(have_sed), "have_anc": bool(have_anc),
                "n_mcparticles": int(n), "schema_version": T.SCHEMA_VERSION}
        meta.update(meta_extra)
        for k in ("DUNESW_VERSION", "SBNDCODE_VERSION", "LARSOFT_VERSION"):
            if os.environ.get(k): meta[k.lower()] = os.environ[k]
        T.write_truth_npz(out, nu, pf, sed_ids=sed_ids, sed_trackq=sed_trackq, anc_map=anc, gen_labels=gen_labels, meta=meta, sed_energy=sed_energy)
        nnu = int((nu[:, T.nu_col("nu_idx")] >= 0).sum()) if len(nu) else 0
        line = ("entry %d r%d s%d e%d: %d nu rows (%d neutrinos), %d particles (%d labelled), %d deposits, %d trackids with charge -> %s"
                % (entry, run, sub, evt, len(nu), nnu, n, int((pf[:, T.pf_col("nu_row")] >= 0).sum()), len(sed_ids), len(sed_trackq), out))
        print(line)
        with open(out[:-4] + ".txt", "w") as f:
            f.write(line + "\n" + json.dumps(meta, sort_keys=True) + "\n")
        ev.next()


if __name__ == "__main__":
    main()
