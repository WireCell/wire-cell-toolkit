#ifndef WIRECELLSIGPROC_PEAKFINDING
#define WIRECELLSIGPROC_PEAKFINDING

#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Logging.h"

namespace WireCell {
    namespace SigProc {

        class PeakFinding {
           public:
            PeakFinding(int fMaxPeaks = 200, double sigma = 1, double threshold = 0.05, bool backgroundRemove = false,
                        int deconIterations = 3, bool markov = true, int averWindow = 3);
            ~PeakFinding();

            int find_peak(Waveform::realseq_t& signal);

            void Clear();

            int GetNPeaks() { return npeaks; };
            double* GetPositionX() { return fPositionX; };
            double* GetPositionY() { return fPositionY; };

           private:
            int fMaxPeaks;
            double sigma;
            double threshold;
            bool backgroundRemove;
            int deconIterations;
            bool markov;
            int averWindow;

            // data ...
            double* source = nullptr;
            int ssize = 0;

            double* destVector = nullptr;
            double* fPositionX = nullptr;
            double* fPositionY = nullptr;

            int npeaks;

            // actual search function ...
            int SearchHighRes();

            Log::logptr_t log;
        };
    }  // namespace SigProc
}  // namespace WireCell
#endif

// Local Variables:
// mode: c++
// c-basic-offset: 2
// End:
