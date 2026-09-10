// doc pdvd/56 T4 -- PR::Fit::reg_flag_u/v/w and TrackFitting's carrier vectors.
//
// Writer-only: no verdict path reads these.  dQ_dx_fit (the single-segment
// path CheckSTM_Michel/TaggerCheckSTM actually call) populates them from the
// dead-channel booleans its own regulariser already computes; dQ_dx_multi_fit
// is untouched (STM never calls it).  There is no existing unit harness that
// exercises dQ_dx_fit itself (it needs a real fitted segment); these cases
// pin the additive struct/member contract the real-arm measurement (doc
// pdvd/60) then relies on.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/PRCommon.h"

using namespace WireCell;
using namespace WireCell::Clus;

TEST_CASE("doc pdvd/56 T4: PR::Fit defaults reg_flag_u/v/w false")
{
    PR::Fit fit;
    CHECK(fit.reg_flag_u == false);
    CHECK(fit.reg_flag_v == false);
    CHECK(fit.reg_flag_w == false);
    // Existing fields untouched by the addition.
    CHECK(fit.dQ == -1);
    CHECK(fit.index == -1);
    CHECK(fit.flag_fix == false);
}

TEST_CASE("doc pdvd/56 T4: PR::Fit reg_flag_u/v/w are independently settable")
{
    PR::Fit fit;
    fit.reg_flag_u = true;
    CHECK(fit.reg_flag_u == true);
    CHECK(fit.reg_flag_v == false);
    CHECK(fit.reg_flag_w == false);
}
