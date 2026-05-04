#include "PeakFinding.h"
#include <iostream>
#include <math.h>

using namespace WireCell;
using namespace WireCell::SigProc;
using namespace std;

PeakFinding::PeakFinding(int fMaxPeaks, double sigma, double threshold, bool backgroundRemove, int deconIterations,
                         bool markov, int averWindow)
  : fMaxPeaks(fMaxPeaks)
  , sigma(sigma)
  , threshold(threshold)
  , backgroundRemove(backgroundRemove)
  , deconIterations(deconIterations)
  , markov(markov)
  , averWindow(averWindow)
  , log(Log::logger("sigproc"))
{
}

PeakFinding::~PeakFinding() { Clear(); }

void PeakFinding::Clear()
{
    delete[] source;
    source = nullptr;
    delete[] destVector;
    destVector = nullptr;
    delete[] fPositionX;
    fPositionX = nullptr;
    delete[] fPositionY;
    fPositionY = nullptr;
}

int PeakFinding::find_peak(Waveform::realseq_t& signal)
{
    Clear();  // free any previous allocations to prevent memory leak
    ssize = int(signal.size());
    source = new double[ssize];
    for (size_t i = 0; i != signal.size(); i++) {
        *(source + i) = signal.at(i);
    }

    destVector = new double[ssize];
    fPositionX = new double[ssize];
    fPositionY = new double[ssize];

    if (ssize >= 2 * sigma + 1) {  // limit the size ...
        npeaks = SearchHighRes();
    }
    else {
        npeaks = 0;
    }

    // fill fPositionY
    for (int i = 0; i < npeaks; i++) {
        int peak_pos1 = std::round(fPositionX[i]);
        fPositionY[i] = source[peak_pos1];
    }

    return npeaks;
}

int PeakFinding::SearchHighRes()
{
    const int PEAK_WINDOW = 1024;

    int i = 0, j = 0, numberIterations = (int) (7 * sigma + 0.5);
    double a = 0, b = 0, c = 0;
    int k = 0, lindex = 0, posit = 0, imin = 0, imax = 0, jmin = 0, jmax = 0, lh_gold = 0, priz = 0;
    double lda = 0, ldb = 0, ldc = 0, area = 0, maximum = 0, maximum_decon = 0;
    int xmin = 0, xmax = 0, l = 0, peak_index = 0, w = 0;
    int size_ext = ssize + 2 * numberIterations, shift = numberIterations, bw = 2;
    double maxch = 0;
    double nom = 0, nip = 0, nim = 0, sp = 0, sm = 0, plocha = 0;
    double m0low = 0, m1low = 0, m2low = 0, l0low = 0, l1low = 0, detlow = 0, av = 0, men = 0;
    if (sigma < 1) {
        log->error("SearchHighRes: invalid sigma {}, must be greater than or equal to 1", sigma);
        return 0;
    }

    if (threshold <= 0 || threshold >= 100) {
        log->error("SearchHighRes: invalid threshold {}, must be positive and less than 100", threshold);
        return 0;
    }

    j = (int) (5.0 * sigma + 0.5);
    if (j >= PEAK_WINDOW / 2) {
        log->error("SearchHighRes: too large sigma: j={}", j);
        return 0;
    }

    if (markov == true) {
        if (averWindow <= 0) {
            log->error("SearchHighRes: averaging window must be positive, got {}", averWindow);
            return 0;
        }
    }

    if (backgroundRemove == true) {
        if (ssize < 2 * numberIterations + 1) {
            log->error("SearchHighRes: too large clipping window: {}", ssize);
            return 0;
        }
    }

    k = int(2 * sigma + 0.5);
    if (k >= 2) {
        for (i = 0; i < k; i++) {
            a = i, b = source[i];
            m0low += 1, m1low += a, m2low += a * a, l0low += b, l1low += a * b;
        }
        detlow = m0low * m2low - m1low * m1low;
        if (detlow != 0)
            l1low = (-l0low * m1low + l1low * m0low) / detlow;

        else
            l1low = 0;
        if (l1low > 0) l1low = 0;
    }

    else {
        l1low = 0;
    }

    i = (int) (7 * sigma + 0.5);
    i = 2 * i;
    double* working_space = new double[7 * (ssize + i)];
    for (j = 0; j < 7 * (ssize + i); j++) working_space[j] = 0;
    for (i = 0; i < size_ext; i++) {
        if (i < shift) {
            a = i - shift;
            working_space[i + size_ext] = source[0] + l1low * a;
            if (working_space[i + size_ext] < 0) working_space[i + size_ext] = 0;
        }

        else if (i >= ssize + shift) {
            a = i - (ssize - 1 + shift);
            working_space[i + size_ext] = source[ssize - 1];
            if (working_space[i + size_ext] < 0) working_space[i + size_ext] = 0;
        }

        else
            working_space[i + size_ext] = source[i - shift];
    }

    if (backgroundRemove == true) {
        for (i = 1; i <= numberIterations; i++) {
            for (j = i; j < size_ext - i; j++) {
                if (markov == false) {
                    a = working_space[size_ext + j];
                    b = (working_space[size_ext + j - i] + working_space[size_ext + j + i]) / 2.0;
                    if (b < a) a = b;

                    working_space[j] = a;
                }

                else {
                    a = working_space[size_ext + j];
                    av = 0;
                    men = 0;
                    for (w = j - bw; w <= j + bw; w++) {
                        if (w >= 0 && w < size_ext) {
                            av += working_space[size_ext + w];
                            men += 1;
                        }
                    }
                    av = av / men;
                    b = 0;
                    men = 0;
                    for (w = j - i - bw; w <= j - i + bw; w++) {
                        if (w >= 0 && w < size_ext) {
                            b += working_space[size_ext + w];
                            men += 1;
                        }
                    }
                    b = b / men;
                    c = 0;
                    men = 0;
                    for (w = j + i - bw; w <= j + i + bw; w++) {
                        if (w >= 0 && w < size_ext) {
                            c += working_space[size_ext + w];
                            men += 1;
                        }
                    }
                    c = c / men;
                    b = (b + c) / 2;
                    if (b < a) av = b;
                    working_space[j] = av;
                }
            }
            for (j = i; j < size_ext - i; j++) working_space[size_ext + j] = working_space[j];
        }
        for (j = 0; j < size_ext; j++) {
            if (j < shift) {
                a = j - shift;
                b = source[0] + l1low * a;
                if (b < 0) b = 0;
                working_space[size_ext + j] = b - working_space[size_ext + j];
            }

            else if (j >= ssize + shift) {
                a = j - (ssize - 1 + shift);
                b = source[ssize - 1];
                if (b < 0) b = 0;
                working_space[size_ext + j] = b - working_space[size_ext + j];
            }

            else {
                working_space[size_ext + j] = source[j - shift] - working_space[size_ext + j];
            }
        }
        for (j = 0; j < size_ext; j++) {
            if (working_space[size_ext + j] < 0) working_space[size_ext + j] = 0;
        }
    }

    for (i = 0; i < size_ext; i++) {
        working_space[i + 6 * size_ext] = working_space[i + size_ext];
    }

    if (markov == true) {
        for (j = 0; j < size_ext; j++) working_space[2 * size_ext + j] = working_space[size_ext + j];
        xmin = 0, xmax = size_ext - 1;
        for (i = 0, maxch = 0; i < size_ext; i++) {
            working_space[i] = 0;
            if (maxch < working_space[2 * size_ext + i]) maxch = working_space[2 * size_ext + i];
            plocha += working_space[2 * size_ext + i];
        }
        if (maxch == 0) {
            delete[] working_space;
            return 0;
        }

        nom = 1;
        working_space[xmin] = 1;
        for (i = xmin; i < xmax; i++) {
            nip = working_space[2 * size_ext + i] / maxch;
            nim = working_space[2 * size_ext + i + 1] / maxch;
            sp = 0, sm = 0;
            for (l = 1; l <= averWindow; l++) {
                if ((i + l) > xmax)
                    a = working_space[2 * size_ext + xmax] / maxch;

                else
                    a = working_space[2 * size_ext + i + l] / maxch;

                b = a - nip;
                if (a + nip <= 0)
                    a = 1;

                else
                    a = sqrt(a + nip);

                b = b / a;
                b = exp(b);
                sp = sp + b;
                if ((i - l + 1) < xmin)
                    a = working_space[2 * size_ext + xmin] / maxch;

                else
                    a = working_space[2 * size_ext + i - l + 1] / maxch;

                b = a - nim;
                if (a + nim <= 0)
                    a = 1;

                else
                    a = sqrt(a + nim);

                b = b / a;
                b = exp(b);
                sm = sm + b;
            }
            a = sp / sm;
            a = working_space[i + 1] = working_space[i] * a;
            nom = nom + a;
        }
        for (i = xmin; i <= xmax; i++) {
            working_space[i] = working_space[i] / nom;
        }
        for (j = 0; j < size_ext; j++) working_space[size_ext + j] = working_space[j] * plocha;
        for (j = 0; j < size_ext; j++) {
            working_space[2 * size_ext + j] = working_space[size_ext + j];
        }
        if (backgroundRemove == true) {
            for (i = 1; i <= numberIterations; i++) {
                for (j = i; j < size_ext - i; j++) {
                    a = working_space[size_ext + j];
                    b = (working_space[size_ext + j - i] + working_space[size_ext + j + i]) / 2.0;
                    if (b < a) a = b;
                    working_space[j] = a;
                }
                for (j = i; j < size_ext - i; j++) working_space[size_ext + j] = working_space[j];
            }
            for (j = 0; j < size_ext; j++) {
                working_space[size_ext + j] = working_space[2 * size_ext + j] - working_space[size_ext + j];
            }
        }
    }
    // deconvolution starts
    area = 0;
    lh_gold = -1;
    posit = 0;
    maximum = 0;
    // generate response vector
    for (i = 0; i < size_ext; i++) {
        lda = (double) i - 3 * sigma;
        lda = lda * lda / (2 * sigma * sigma);
        j = (int) (1000 * exp(-lda));
        lda = j;
        if (lda != 0) lh_gold = i + 1;

        working_space[i] = lda;
        area = area + lda;
        if (lda > maximum) {
            maximum = lda;
            posit = i;
        }
    }
    // read source vector
    for (i = 0; i < size_ext; i++) working_space[2 * size_ext + i] = fabs(working_space[size_ext + i]);
    // create matrix at*a(vector b)
    i = lh_gold - 1;
    if (i > size_ext) i = size_ext;

    imin = -i, imax = i;
    for (i = imin; i <= imax; i++) {
        lda = 0;
        jmin = 0;
        if (i < 0) jmin = -i;
        jmax = lh_gold - 1 - i;
        if (jmax > (lh_gold - 1)) jmax = lh_gold - 1;

        for (j = jmin; j <= jmax; j++) {
            ldb = working_space[j];
            ldc = working_space[i + j];
            lda = lda + ldb * ldc;
        }
        working_space[size_ext + i - imin] = lda;
    }
    // create vector p
    i = lh_gold - 1;
    imin = -i, imax = size_ext + i - 1;
    for (i = imin; i <= imax; i++) {
        lda = 0;
        for (j = 0; j <= (lh_gold - 1); j++) {
            ldb = working_space[j];
            k = i + j;
            if (k >= 0 && k < size_ext) {
                ldc = working_space[2 * size_ext + k];
                lda = lda + ldb * ldc;
            }
        }
        working_space[4 * size_ext + i - imin] = lda;
    }
    // move vector p
    for (i = imin; i <= imax; i++) working_space[2 * size_ext + i - imin] = working_space[4 * size_ext + i - imin];
    // initialization of resulting vector
    for (i = 0; i < size_ext; i++) working_space[i] = 1;
    // START OF ITERATIONS
    for (lindex = 0; lindex < deconIterations; lindex++) {
        for (i = 0; i < size_ext; i++) {
            if (fabs(working_space[2 * size_ext + i]) > 0.00001 && fabs(working_space[i]) > 0.00001) {
                lda = 0;
                jmin = lh_gold - 1;
                if (jmin > i) jmin = i;

                jmin = -jmin;
                jmax = lh_gold - 1;
                if (jmax > (size_ext - 1 - i)) jmax = size_ext - 1 - i;

                for (j = jmin; j <= jmax; j++) {
                    ldb = working_space[j + lh_gold - 1 + size_ext];
                    ldc = working_space[i + j];
                    lda = lda + ldb * ldc;
                }
                ldb = working_space[2 * size_ext + i];
                if (lda != 0)
                    lda = ldb / lda;

                else
                    lda = 0;

                ldb = working_space[i];
                lda = lda * ldb;
                working_space[3 * size_ext + i] = lda;
            }
        }
        for (i = 0; i < size_ext; i++) {
            working_space[i] = working_space[3 * size_ext + i];
        }
    }
    // shift resulting spectrum
    for (i = 0; i < size_ext; i++) {
        lda = working_space[i];
        j = i + posit;
        j = j % size_ext;
        working_space[size_ext + j] = lda;
    }
    // write back resulting spectrum
    maximum = 0, maximum_decon = 0;
    j = lh_gold - 1;
    for (i = 0; i < size_ext - j; i++) {
        if (i >= shift && i < ssize + shift) {
            working_space[i] = area * working_space[size_ext + i + j];
            if (maximum_decon < working_space[i]) maximum_decon = working_space[i];
            if (maximum < working_space[6 * size_ext + i]) maximum = working_space[6 * size_ext + i];
        }

        else
            working_space[i] = 0;
    }
    lda = 1;
    if (lda > threshold) lda = threshold;
    lda = lda / 100;

    // searching for peaks in deconvolved spectrum
    for (i = 1; i < size_ext - 1; i++) {
        if (working_space[i] > working_space[i - 1] && working_space[i] > working_space[i + 1]) {
            if (i >= shift && i < ssize + shift) {
                if (working_space[i] > lda * maximum_decon &&
                    working_space[6 * size_ext + i] > threshold * maximum / 100.0) {
                    for (j = i - 1, a = 0, b = 0; j <= i + 1; j++) {
                        a += (double) (j - shift) * working_space[j];
                        b += working_space[j];
                    }
                    a = a / b;
                    if (a < 0) a = 0;

                    if (a >= ssize) a = ssize - 1;
                    if (peak_index == 0) {
                        fPositionX[0] = a;
                        peak_index = 1;
                    }

                    else {
                        for (j = 0, priz = 0; j < peak_index && priz == 0; j++) {
                            if (working_space[6 * size_ext + shift + (int) a] >
                                working_space[6 * size_ext + shift + (int) fPositionX[j]])
                                priz = 1;
                        }
                        if (priz == 0) {
                            if (j < fMaxPeaks) {
                                fPositionX[j] = a;
                            }
                        }

                        else {
                            for (k = peak_index; k >= j; k--) {
                                if (k < fMaxPeaks) {
                                    fPositionX[k] = fPositionX[k - 1];
                                }
                            }
                            fPositionX[j - 1] = a;
                        }
                        if (peak_index < fMaxPeaks) peak_index += 1;
                    }
                }
            }
        }
    }

    for (i = 0; i < ssize; i++) destVector[i] = working_space[i + shift];
    delete[] working_space;
    // fNPeaks = peak_index;
    // if(peak_index == fMaxPeaks){
    //   log->debug("SearchHighRes: peak buffer full with {} peaks, nticks: {}, theshold: {}",
    //              fMaxPeaks, ssize, threshold);

    // }
    return peak_index;
}

// Local Variables:
// mode: c++
// c-basic-offset: 2
// End:
