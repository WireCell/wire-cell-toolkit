// This is used to make a figure in the dnn-roi paper.
// Please do not change.

#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/RaySolving.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Logging.h"

#include "TCanvas.h"
#include "TMarker.h"
#include "TText.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TArrow.h"
#include "TH1F.h"

#include "WireCellUtil/Graph.h"

#include <math.h>

#include <random>

using namespace WireCell;
using namespace WireCell::Waveform;
using namespace WireCell::RayGrid;
using namespace std;
using spdlog::info;
using spdlog::warn;

const int ndepos = 10;
const int neles = 10;
const double pitch_magnitude = 5;
const double gaussian = 3;
const double border = 10;
const double width = 100;
const double height = 100;

// local helper codes
#include "raygrid.h"
#include "raygrid_draw.h"

int main(int argc, char* argv[])
{
    Printer print(argv[0]);

    auto raypairs = make_raypairs(width, height, pitch_magnitude);
    //const int nlayers = raypairs.size();

    Coordinates coords(raypairs);

    // This is used to make a figure in the dnn-roi paper.
    // Please do not change.
    auto frame = print.canvas.DrawFrame(-110, -25, 75, 125);
    draw_raygrid(print, coords, raypairs, frame);
    return 0;
}
