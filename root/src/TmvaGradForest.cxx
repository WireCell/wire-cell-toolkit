#include "WireCellRoot/TmvaGradForest.h"
#include "WireCellUtil/Exceptions.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string_view>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace WireCell;
using WireCell::Root::TmvaGradForest;

namespace {

    using sv = std::string_view;

    [[noreturn]] void bail(const std::string& what)
    {
        THROW(ValueError() << errmsg{"TmvaGradForest: " + what});
    }

    // The attribute list of one start tag: everything between the tag name
    // and the closing '>' (or "/>").  TMVA/xgboost2TMVA never put '>' inside
    // an attribute value, and neither does any XML writer that escapes.
    struct Tag {
        sv attrs;
        bool selfclose{false};
        size_t end{0};  // offset just past the '>'
    };

    Tag read_tag(sv doc, size_t at)
    {
        // `at` points at '<'.  Skip the name.
        size_t p = at + 1;
        while (p < doc.size() && doc[p] != ' ' && doc[p] != '>' && doc[p] != '/') ++p;
        size_t a0 = p;
        size_t gt = doc.find('>', p);
        if (gt == sv::npos) bail("unterminated tag");
        Tag t;
        t.selfclose = gt > 0 && doc[gt - 1] == '/';
        t.attrs = doc.substr(a0, (t.selfclose ? gt - 1 : gt) - a0);
        t.end = gt + 1;
        return t;
    }

    // Find name="value" inside an attribute list; empty sv + false if absent.
    // Matches whole attribute names only (so "Cut" never matches "NCut...").
    bool attr(sv attrs, sv name, sv& out)
    {
        size_t p = 0;
        while (p < attrs.size()) {
            // skip whitespace
            while (p < attrs.size() && (attrs[p] == ' ' || attrs[p] == '\n' || attrs[p] == '\t' || attrs[p] == '\r')) ++p;
            if (p >= attrs.size()) break;
            size_t eq = attrs.find('=', p);
            if (eq == sv::npos) break;
            sv key = attrs.substr(p, eq - p);
            size_t q1 = attrs.find('"', eq);
            if (q1 == sv::npos) break;
            size_t q2 = attrs.find('"', q1 + 1);
            if (q2 == sv::npos) break;
            if (key == name) {
                out = attrs.substr(q1 + 1, q2 - q1 - 1);
                return true;
            }
            p = q2 + 1;
        }
        return false;
    }

    sv need(sv attrs, sv name, const char* where)
    {
        sv v;
        if (!attr(attrs, name, v)) bail(std::string("missing attribute '") + std::string(name) + "' on " + where);
        return v;
    }

    // Tools::ReadAttr(float&): value = atof(val).  atof == strtod(val, 0);
    // the attribute value is followed by '"' in the document, at which
    // strtod stops, so it can be called in place without a copy.
    float read_float(sv v) { return static_cast<float>(std::strtod(v.data(), nullptr)); }

    // Tools::ReadAttr(int&/short&): atoi(val) == (int) strtol(val, 0, 10).
    int read_int(sv v) { return static_cast<int>(std::strtol(v.data(), nullptr, 10)); }

    // Tools::ReadAttr<Bool_t>: std::stringstream(val) >> value, i.e. the
    // integer form ("0"/"1"); on a failed extraction the stream leaves the
    // Bool_t as it was.  DecisionTreeNode's fCutType defaults to kTRUE.
    bool read_bool(sv v)
    {
        std::stringstream s{std::string(v)};
        bool value = true;
        s >> value;
        return value;
    }

    // Text between the start tag at `at` and the next '<' (element content).
    sv content_after(sv doc, size_t tag_end)
    {
        size_t lt = doc.find('<', tag_end);
        if (lt == sv::npos) return {};
        return doc.substr(tag_end, lt - tag_end);
    }

    // Case-insensitive, whitespace-trimmed compare (TMVA's option values are
    // matched after ToLower()/trim in ReadOptionsFromXML).
    std::string norm(sv v)
    {
        size_t b = 0, e = v.size();
        while (b < e && std::isspace(static_cast<unsigned char>(v[b]))) ++b;
        while (e > b && std::isspace(static_cast<unsigned char>(v[e - 1]))) --e;
        std::string s(v.substr(b, e - b));
        for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    }

}  // namespace

void TmvaGradForest::load(const std::string& path, const std::vector<std::string>& names)
{
    // Map the file read-only instead of copying it into a string: the
    // nue combiner is 199 MB, and the copy would be the largest transient
    // of the whole load.  The parser only needs a string_view.
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) bail("cannot open " + path);
    struct stat st;
    if (::fstat(fd, &st) != 0 || st.st_size <= 0) {
        ::close(fd);
        bail("empty or unreadable file " + path);
    }
    const size_t size = static_cast<size_t>(st.st_size);
    void* map = ::mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0);
    ::close(fd);
    if (map == MAP_FAILED) bail("mmap failed on " + path);
    try {
        parse(sv(static_cast<const char*>(map), size), names);
    }
    catch (...) {
        ::munmap(map, size);
        throw;
    }
    ::munmap(map, size);
}

void TmvaGradForest::load_string(const std::string& xml, const std::vector<std::string>& names)
{
    parse(sv(xml), names);
}

void TmvaGradForest::parse(sv doc, const std::vector<std::string>& names)
{
    m_nodes.clear();
    m_roots.clear();
    m_nvars = 0;

    // ---- <MethodSetup Method="BDT::BDT"> --------------------------------
    size_t ms = doc.find("<MethodSetup");
    if (ms == sv::npos) bail("no <MethodSetup> element");
    {
        Tag t = read_tag(doc, ms);
        sv method = need(t.attrs, "Method", "<MethodSetup>");
        if (method != "BDT::BDT") bail("unsupported Method '" + std::string(method) + "' (need BDT::BDT)");
    }

    // ---- <Variables NVar=..> in caller order (ReadVariablesFromXML) ------
    size_t vs = doc.find("<Variables", ms);
    if (vs == sv::npos) bail("no <Variables> element");
    {
        Tag t = read_tag(doc, vs);
        const int nvar = read_int(need(t.attrs, "NVar", "<Variables>"));
        if (nvar < 0 || static_cast<size_t>(nvar) != names.size()) {
            bail("variable count mismatch: file has " + std::to_string(nvar) + ", caller registered " +
                 std::to_string(names.size()));
        }
        size_t vend = doc.find("</Variables>", t.end);
        if (vend == sv::npos) bail("unterminated <Variables>");
        size_t p = t.end;
        for (int i = 0; i < nvar; ++i) {
            size_t v = doc.find("<Variable ", p);
            if (v == sv::npos || v > vend) bail("fewer <Variable> elements than NVar");
            Tag vt = read_tag(doc, v);
            sv expr = need(vt.attrs, "Expression", "<Variable>");
            if (expr != sv(names[i])) {
                bail("variable " + std::to_string(i) + " mismatch: file '" + std::string(expr) + "', caller '" +
                     names[i] + "'");
            }
            p = vt.end;
        }
        m_nvars = static_cast<size_t>(nvar);
    }

    // ---- refuse anything outside the mirrored subset ----------------------
    {
        size_t tr = doc.find("<Transformations", ms);
        if (tr != sv::npos) {
            Tag t = read_tag(doc, tr);
            sv n;
            if (attr(t.attrs, "NTransformations", n) && read_int(n) != 0) {
                bail("input-variable transformations present (NTransformations=" + std::string(n) + ")");
            }
        }
        // Options: BoostType must be Grad (MethodBDT::PrivateGetMvaValue takes
        // the GradBoost branch only then); DoPreselection must not be on.
        std::string boost;
        size_t op = doc.find("<Options", ms);
        if (op != sv::npos) {
            size_t oend = doc.find("</Options>", op);
            if (oend == sv::npos) bail("unterminated <Options>");
            size_t p = op;
            while (true) {
                size_t o = doc.find("<Option ", p);
                if (o == sv::npos || o > oend) break;
                Tag ot = read_tag(doc, o);
                sv name = need(ot.attrs, "name", "<Option>");
                std::string val = norm(content_after(doc, ot.end));
                if (name == "BoostType") boost = val;
                if (name == "DoPreselection" && val != "false" && val != "0") bail("DoPreselection is on");
                p = ot.end;
            }
        }
        if (boost != "grad") bail("BoostType '" + boost + "' is not Grad");
    }

    // ---- <Weights NTrees=.. AnalysisType="1"> ----------------------------
    size_t ws = doc.find("<Weights", ms);
    if (ws == sv::npos) bail("no <Weights> element");
    int ntrees = 0;
    {
        Tag t = read_tag(doc, ws);
        sv dummy;
        if (attr(t.attrs, "PreselectionLowBkgVar0", dummy)) bail("preselection cuts present");
        ntrees = read_int(need(t.attrs, "NTrees", "<Weights>"));
        sv at;
        if (!attr(t.attrs, "AnalysisType", at) && !attr(t.attrs, "TreeType", at)) {
            bail("<Weights> has neither AnalysisType nor TreeType");
        }
        // Types::kRegression == 1: DecisionTree::CheckEvent returns the leaf
        // response.  Any other tree type returns nType/purity -- not mirrored.
        if (read_int(at) != 1) bail("Weights AnalysisType " + std::string(at) + " is not 1 (regression leaves)");
        ws = t.end;
    }
    size_t wend = doc.find("</Weights>", ws);
    if (wend == sv::npos) bail("unterminated <Weights>");

    // ---- trees -------------------------------------------------------------
    // MethodBDT::ReadWeightsFromXML walks every child of <Weights> as a tree;
    // BinaryTree::ReadXML reads the tree's first child as the root node;
    // Node::ReadXML recurses, attaching children by their pos attribute.
    m_roots.reserve(static_cast<size_t>(ntrees));
    std::vector<int32_t> stack;  // open <Node> elements
    size_t p = ws;
    int tree_root = -1;  // index of the current tree's root, -1 between trees
    bool in_tree = false;
    while (true) {
        size_t lt = doc.find('<', p);
        if (lt == sv::npos || lt >= wend) break;
        sv rest = doc.substr(lt);
        if (rest.compare(0, 12, "<BinaryTree ") == 0 || rest.compare(0, 12, "<BinaryTree>") == 0) {
            Tag t = read_tag(doc, lt);
            if (in_tree) bail("nested <BinaryTree>");
            in_tree = true;
            tree_root = -1;
            stack.clear();
            p = t.end;
            continue;
        }
        if (rest.compare(0, 13, "</BinaryTree>") == 0) {
            if (!in_tree) bail("stray </BinaryTree>");
            if (!stack.empty()) bail("unclosed <Node> at </BinaryTree>");
            if (tree_root < 0) bail("<BinaryTree> without a root node");
            m_roots.push_back(tree_root);
            in_tree = false;
            p = lt + 13;
            continue;
        }
        if (rest.compare(0, 6, "<Node ") == 0) {
            if (!in_tree) bail("<Node> outside <BinaryTree>");
            Tag t = read_tag(doc, lt);
            Node n;
            sv v;
            if (attr(t.attrs, "NCoef", v) && read_int(v) != 0) bail("Fisher-cut node (NCoef != 0)");
            n.ivar = read_int(need(t.attrs, "IVar", "<Node>"));
            n.cut = read_float(need(t.attrs, "Cut", "<Node>"));
            n.ctype = read_bool(need(t.attrs, "cType", "<Node>"));
            n.res = attr(t.attrs, "res", v) ? read_float(v) : -99.f;  // DecisionTreeNode default fResponse
            n.ntype = read_int(need(t.attrs, "nType", "<Node>"));
            sv pos = need(t.attrs, "pos", "<Node>");
            if (n.ntype == 0 && (n.ivar < 0 || static_cast<size_t>(n.ivar) >= m_nvars)) {
                bail("intermediate node selects variable " + std::to_string(n.ivar) + " of " + std::to_string(m_nvars));
            }
            const int32_t idx = static_cast<int32_t>(m_nodes.size());
            m_nodes.push_back(n);
            if (stack.empty()) {
                if (tree_root >= 0) bail("two root nodes in one <BinaryTree>");
                tree_root = idx;
            }
            else {
                Node& parent = m_nodes[stack.back()];
                if (pos == "l") {
                    if (parent.left >= 0) bail("node with two left children");
                    parent.left = idx;
                }
                else if (pos == "r") {
                    if (parent.right >= 0) bail("node with two right children");
                    parent.right = idx;
                }
                else {
                    bail("child node with pos '" + std::string(pos) + "' (neither left nor right)");
                }
            }
            if (!t.selfclose) stack.push_back(idx);
            p = t.end;
            continue;
        }
        if (rest.compare(0, 7, "</Node>") == 0) {
            if (stack.empty()) bail("stray </Node>");
            stack.pop_back();
            p = lt + 7;
            continue;
        }
        bail("unexpected markup inside <Weights>: " + std::string(rest.substr(0, 20)));
    }
    if (in_tree) bail("unterminated <BinaryTree>");
    if (static_cast<int>(m_roots.size()) != ntrees) {
        bail("NTrees=" + std::to_string(ntrees) + " but " + std::to_string(m_roots.size()) + " trees found");
    }
    // Every intermediate node must have both children (DT::CheckEvent would
    // hit a null pointer and abort otherwise).
    for (const auto& n : m_nodes) {
        if (n.ntype == 0 && (n.left < 0 || n.right < 0)) bail("intermediate node missing a child");
    }
}

double TmvaGradForest::evaluate(const float* values) const
{
    // Reader::EvaluateMVA: NaN in any input variable => -999.
    for (size_t i = 0; i < m_nvars; ++i) {
        if (std::isnan(values[i])) return -999;
    }
    // MethodBDT::GetGradBoostMVA over DecisionTree::CheckEvent(ev, kFALSE)
    // on regression trees (leaf response), nodes stepped by
    // DecisionTreeNode::GoesRight.
    double sum = 0;
    for (const int32_t root : m_roots) {
        const Node* cur = &m_nodes[root];
        while (cur->ntype == 0) {
            bool result = (values[cur->ivar] >= cur->cut);
            if (!cur->ctype) result = !result;
            cur = &m_nodes[result ? cur->right : cur->left];
        }
        sum += cur->res;
    }
    return 2.0 / (1.0 + exp(-2.0 * sum)) - 1;
}
