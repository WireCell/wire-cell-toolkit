/** Simple implementations of the trigger object sets. */

#ifndef WIRECELLAUX_SIMPLETRIGGERSETS
#define WIRECELLAUX_SIMPLETRIGGERSETS

#include "WireCellIface/ITriggerCandidateSet.h"

namespace WireCell::Aux {

    class SimpleTriggerPrimitiveSet : public ITriggerPrimitiveSet {
       public:
        SimpleTriggerPrimitiveSet(int ident, int detector, double time_begin, double time_end,
                                  primitives_t primitives = {})
          : m_ident(ident), m_detector(detector), m_begin(time_begin), m_end(time_end)
          , m_objects(std::move(primitives))
        {
        }
        virtual ~SimpleTriggerPrimitiveSet();
        virtual int ident() const { return m_ident; }
        virtual int detector() const { return m_detector; }
        virtual double time_begin() const { return m_begin; }
        virtual double time_end() const { return m_end; }
        virtual const primitives_t& primitives() const { return m_objects; }

       private:
        int m_ident, m_detector;
        double m_begin, m_end;
        primitives_t m_objects;
    };

    class SimpleTriggerActivitySet : public ITriggerActivitySet {
       public:
        SimpleTriggerActivitySet(int ident, int detector, double time_begin, double time_end,
                                 activities_t activities = {})
          : m_ident(ident), m_detector(detector), m_begin(time_begin), m_end(time_end)
          , m_objects(std::move(activities))
        {
        }
        virtual ~SimpleTriggerActivitySet();
        virtual int ident() const { return m_ident; }
        virtual int detector() const { return m_detector; }
        virtual double time_begin() const { return m_begin; }
        virtual double time_end() const { return m_end; }
        virtual const activities_t& activities() const { return m_objects; }

       private:
        int m_ident, m_detector;
        double m_begin, m_end;
        activities_t m_objects;
    };

    class SimpleTriggerCandidateSet : public ITriggerCandidateSet {
       public:
        SimpleTriggerCandidateSet(int ident, int detector, double time_begin, double time_end,
                                  candidates_t candidates = {})
          : m_ident(ident), m_detector(detector), m_begin(time_begin), m_end(time_end)
          , m_objects(std::move(candidates))
        {
        }
        virtual ~SimpleTriggerCandidateSet();
        virtual int ident() const { return m_ident; }
        virtual int detector() const { return m_detector; }
        virtual double time_begin() const { return m_begin; }
        virtual double time_end() const { return m_end; }
        virtual const candidates_t& candidates() const { return m_objects; }

       private:
        int m_ident, m_detector;
        double m_begin, m_end;
        candidates_t m_objects;
    };

}  // namespace WireCell::Aux

#endif
