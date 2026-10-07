// ai-helper issue 35: the record of one DL-vertex network call (PRDlVtxDump.h).
// Pin the contract the writer and the standalone replay rely on: a fresh record
// is a production-pass, status-0 call with no decision, no payload and every row
// index at -1 (-1 = "not a vertex of this cloud", the writer copies it verbatim).
#include "WireCellClus/PRDlVtxDump.h"
#include "WireCellUtil/doctest.h"

using WireCell::Clus::PR::DlVtxCall;

TEST_CASE("dlvtx dump: DlVtxCall defaults")
{
    DlVtxCall c;
    CHECK(c.pass == "prod");
    CHECK(c.status == 0);
    CHECK(c.top_k == 0);
    CHECK_FALSE(c.rerank);
    CHECK(c.x.empty()); CHECK(c.y.empty()); CHECK(c.z.empty()); CHECK(c.q.empty());
    CHECK(c.n_vertex_rows == 0);
    CHECK_FALSE(c.cloud_no_exclusion);
    CHECK(c.payload.empty());
    CHECK_FALSE(c.payload_from_off);
    CHECK(c.n_off_voxels == 0);
    CHECK_FALSE(c.trad_valid);   CHECK(c.trad_row == -1);
    CHECK_FALSE(c.rerank_valid); CHECK(c.rerank_row == -1);
    CHECK_FALSE(c.accepted);     CHECK(c.dl_row == -1);
    CHECK_FALSE(c.dual_transferred);
    CHECK_FALSE(c.two_end_veto);
    CHECK_FALSE(c.hint_valid);
}

TEST_CASE("dlvtx dump: the cloud columns stay parallel")
{
    DlVtxCall c;
    c.x = {1.f, 2.f}; c.y = {3.f, 4.f}; c.z = {5.f, 6.f}; c.q = {7.f, 8.f};
    c.n_vertex_rows = 1;
    CHECK(c.x.size() == c.y.size());
    CHECK(c.x.size() == c.z.size());
    CHECK(c.x.size() == c.q.size());
    CHECK(c.n_vertex_rows <= static_cast<int>(c.x.size()));
}
