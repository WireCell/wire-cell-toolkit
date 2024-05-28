#include <random>

#include <json/json.h>

#include <iostream>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <string>

#include <fstream>
using namespace std;

int main(int argc, char* argv[])
{
    string name = argv[0];

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    auto uni = std::bind(distribution, generator);

    cout << "Generating data\n";
    const int nbins = 500;
    const int narrays = 100;
    double data[narrays][nbins];
    Json::Value jroot;
    for (int ind = 0; ind < narrays; ++ind) {
        Json::Value jarr;
        for (int ibin = 0; ibin < nbins; ++ibin) {
            auto number = uni();
            data[ind][ibin] = number;
            jarr.append(number);
        }
        jroot[ind] = jarr;
    }

    {
        string fname = name + ".bin";
        cout << "writing " << fname << endl;
        fstream myfile(fname.c_str(), ios::out | ios::binary);
        myfile.write((char*) data, nbins * narrays * sizeof(double));
        myfile.close();
    }

    {
        string fname = name + ".json";
        cout << "writing " << fname << endl;
        ofstream jout(fname.c_str());
        jout << jroot;
    }

    {
        string fname = name + ".json.bz2";
        cout << "writing " << fname << endl;
        std::fstream jbz2file(fname.c_str(), std::ios::binary | std::ios::out);
        boost::iostreams::filtering_stream<boost::iostreams::output> outfilt;
        outfilt.push(boost::iostreams::bzip2_compressor());
        outfilt.push(jbz2file);
        outfilt << jroot;
    }

    Json::Value jroot2;
    {
        string fname = name + ".json.bz2";
        cout << "reading " << fname << endl;
        std::fstream jbz2file(fname.c_str(), std::ios::binary | std::ios::in);
        boost::iostreams::filtering_stream<boost::iostreams::input> infilt;
        infilt.push(boost::iostreams::bzip2_decompressor());
        infilt.push(jbz2file);
        infilt >> jroot2;
    }

    {
        string fname = name + "2.json";
        cout << "writing " << fname << endl;
        ofstream jout2(fname.c_str());
        jout2 << jroot2;
    }

    return 0;
}
