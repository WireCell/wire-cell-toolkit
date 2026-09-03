/** doc pdvd/26 round 2: a pass budget for the PR tail's `while (flag_continue)`
 *  examiners.
 *
 *  Every examiner in find_proto_vertex (examine_structure_2/3, examine_vertices,
 *  examine_partial_identical_segments, examine_structure_final_1/2/3,
 *  eliminate_short_vertex_activities, the two shower_clustering_with_nv_*
 *  loops) re-runs until a full pass changes nothing.  Nothing bounds the
 *  number of passes; a state that is a fixed point of one pass -- doc pdvd/26,
 *  PDVD 039349/14 and /53, one vertex cloned in place per 0.42 s pass for
 *  48 min -- looks exactly like a live job.  This counter turns that into a
 *  WARN and a finished event.  The budget is far above the observed maximum
 *  (doc pdvd/26 sec 7 census: the DEBUG line below over 120 PDVD, 67 SBND and
 *  35 uBooNE events), so on every event that terminates today the counter
 *  is a no-op and the outputs are byte-identical.
 *
 *  Usage, at the loop:
 *      ExaminerPassCounter epb("examine_vertices", cluster.get_cluster_id());
 *      while (flag_continue) {
 *          if (epb.exceeded()) break;
 *          flag_continue = false;
 *          ...
 *  The destructor emits one DEBUG line per call ("examiner passes: <who>
 *  cluster <id> passes=<n>") -- the census input.
 */
#pragma once

namespace WireCell::Clus::PR {

    inline constexpr int kExaminerPassBudget = 1000;

    struct ExaminerPassCounter {
        const char* who;
        int cluster_id;
        int passes = 0;
        ExaminerPassCounter(const char* w, int cid) : who(w), cluster_id(cid) {}
        /// Call once per pass, first thing in the while body.  Increments the
        /// pass count; true (and a WARN) once the budget is exceeded.
        bool exceeded();
        ~ExaminerPassCounter();
    };

}  // namespace WireCell::Clus::PR
