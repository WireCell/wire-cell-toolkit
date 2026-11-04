#include "WireCellClus/PerEventName.h"
#include "WireCellUtil/doctest.h"

using WireCell::Clus::per_event_name;
using WireCell::Clus::per_event_templated;

TEST_CASE("per_event_name: the original ident idiom, literals and the RSE tokens")
{
    CHECK(per_event_name("pr_evt%1%", 7, 1, 2, 7) == "pr_evt7");
    CHECK(per_event_name("pr_evt%1%/tracking-pr.root", 7, 1, 2, 7) == "pr_evt7/tracking-pr.root");
    CHECK(per_event_name("tracking-pr.root", 7, 1, 2, 7) == "tracking-pr.root");
    CHECK(per_event_name("pr_evt{run}_{subrun}_{event}", 1, 2100126, 29, 1) == "pr_evt2100126_29_1");
    CHECK(per_event_name("pr_evt{ident}/x%1%", 7, 1, 2, 3) == "pr_evt7/x7");
    CHECK(per_event_name("{event}{event}", 7, 1, 2, 3) == "33");
    CHECK(per_event_templated("pr_evt%1%"));
    CHECK(per_event_templated("pr_evt{event}"));
    CHECK_FALSE(per_event_templated("mabc-pr.zip"));
    // PR 544 review, item 10: only the known tokens make a template.
    CHECK(per_event_templated("pr_evt{run}_{subrun}_{event}"));
    CHECK(per_event_templated("x{ident}.root"));
    CHECK_FALSE(per_event_templated("weird{name}.zip"));
    CHECK(per_event_name("weird{name}.zip", 7, 1, 2, 3) == "weird{name}.zip");
}
