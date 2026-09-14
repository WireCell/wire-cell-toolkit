// sbnd_xin/docs/110 -- PR::reset_shower_id_counter() restarts the process-wide
// shower-id sequence.
//
// Shower ids come from one static counter (PRShower.cxx), so a process that
// streams several events (the SBND group mode) numbered event N's showers after
// event N-1's, and the PrDisplayDump calib json of every event but the first
// disagreed with the one-event-per-process job.  MultiAlgBlobClustering calls the
// reset at each event start when reset_shower_ids_per_event is on.
//
// Making the reset a no-op must fail the first case.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRShower.h"

#include <memory>

using namespace WireCell::Clus;

TEST_CASE("doc110 shower id counter: a reset restarts the sequence at 0")
{
    PR::Graph graph;
    auto a = std::make_shared<PR::Shower>(graph);
    auto b = std::make_shared<PR::Shower>(graph);
    CHECK(b->get_shower_id() == a->get_shower_id() + 1);  // monotonic within a process

    PR::reset_shower_id_counter();
    auto c = std::make_shared<PR::Shower>(graph);
    auto d = std::make_shared<PR::Shower>(graph);
    CHECK(c->get_shower_id() == 0);
    CHECK(d->get_shower_id() == 1);
}

TEST_CASE("doc110 shower id counter: without a reset the next event keeps counting")
{
    PR::Graph graph;
    PR::reset_shower_id_counter();
    {
        auto first_event = std::make_shared<PR::Shower>(graph);
        CHECK(first_event->get_shower_id() == 0);
    }
    // A second "event" in the same process, no reset: the legacy numbering the
    // knob-off path keeps.
    auto second_event = std::make_shared<PR::Shower>(graph);
    CHECK(second_event->get_shower_id() == 1);
}
