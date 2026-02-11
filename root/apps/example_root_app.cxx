// Example ROOT application
// This file demonstrates how to create an application with ROOT dependencies

#include <iostream>
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"

int main(int argc, char* argv[]) {

    std::cout << "Example ROOT Application" << std::endl;

    // Create a simple histogram
    TH1F* hist = new TH1F("hist", "Example Histogram", 100, -5, 5);

    // Fill with random Gaussian data
    hist->FillRandom("gaus", 10000);

    // Save to file
    const char* filename = "example_output.root";
    if (argc > 1) {
        filename = argv[1];
    }

    TFile* outfile = new TFile(filename, "RECREATE");
    hist->Write();
    outfile->Close();

    std::cout << "Saved histogram to " << filename << std::endl;

    delete hist;
    delete outfile;

    return 0;
}