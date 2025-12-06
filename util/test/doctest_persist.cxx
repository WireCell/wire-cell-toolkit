#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/doctest.h"
#include <boost/filesystem/fstream.hpp>
#include <cstdlib>

using spdlog::debug;
using namespace WireCell;

TEST_CASE("persist checks")
{
    CHECK(Persist::exists("/etc"));
    CHECK(Persist::exists("/etc/hosts"));

    setenv("WIRECELL_PATH", "/etc:/usr:/var", 1);
    std::string etchosts = Persist::resolve("hosts");
    // std::cerr << etchosts << std::endl;
    CHECK(etchosts == "/etc/hosts");

    std::string dne = Persist::resolve("this_file-really_should-not_exist");
    CHECK(dne.empty());

}


TEST_CASE("persist tempdir")
{
    boost::filesystem::path path; // persist beyond each scope

    {
        Persist::TempDir p;
        path = p.path;
        debug("tempdir: {}", path.native());
        REQUIRE(boost::filesystem::exists(path));
    }
    REQUIRE(! boost::filesystem::exists(path));

    {
        Persist::TempDir p;
        p.keep = true;
        path = p.path;
        debug("tempdir: {}", path.native());
        REQUIRE(boost::filesystem::exists(path));
    }
    REQUIRE(boost::filesystem::exists(path));
    boost::filesystem::remove_all(path);

    {
        Persist::TempDir p("foo-%%%", true);
        path = p.path;
        debug("tempdir: {}", path.native());
        REQUIRE(boost::filesystem::exists(path));
        auto tdir = boost::filesystem::temp_directory_path();
        REQUIRE(boost::filesystem::exists(tdir / path.filename()));
    }
    REQUIRE(! boost::filesystem::exists(path));

}


TEST_CASE("persist json")
{
    boost::filesystem::path path; // persist beyond each scope
    {
        auto tdir = Persist::TempDir();
        path = tdir.path;

        const int nbins = 500;
        const int narrays = 100;
        Json::Value jroot;
        for (int ind = 0; ind < narrays; ++ind) {
            Json::Value jarr;
            for (int ibin = 0; ibin < nbins; ++ibin) {
                auto number = (ind+1)*(ibin+1);
                jarr.append(number);
            }
            jroot[ind] = jarr;
        }


        // dump non-pretty in plain text and bzip2'ed to see size diff
        Persist::dump((path / "dump.json").native(), jroot);
        REQUIRE(boost::filesystem::exists(path / "dump.json"));
        Persist::dump((path / "dump-ugly.json.bz2").native(), jroot);

        // dump pretty, read back and redump so can do a diff
        Persist::dump((path / "dump-pretty.json").native(), jroot, true);
        Persist::dump((path / "dump-pretty.json.bz2").native(), jroot, true);

        Json::Value jroot2 = Persist::load((path / "dump-pretty.json").native());
        Json::Value jroot3 = Persist::load((path / "dump-pretty.json.bz2").native());

        Persist::dump((path / "dump-pretty2.json").native(), jroot2, true);
        Persist::dump((path / "dump-pretty3.json").native(), jroot3, true);

        auto jroot3b = Persist::load((path / "dump-pretty3.json").native());
        CHECK(jroot3b == jroot3);

        // issue #143 says a WIRECELL_PATH with a non-directory causes a crash
        // so include a file and an empty entry.  
        std::string wpath = (path / "dump-pretty2.json").native() + "::" + path.native();
        setenv("WIRECELL_PATH", wpath.c_str(), 1);
        auto full = Persist::resolve("dump-pretty3.json");
        debug("resolved {}", full);
        REQUIRE(full != "");

        // try a jsonnet parse
        auto psub = path / "sub";
        boost::filesystem::create_directories(psub);

        Persist::dump((psub / "dump.json").native(), jroot);
        auto jpath = path / "something.jsonnet";
        boost::filesystem::ofstream jcode(jpath);
        jcode << "local blah = import \"dump.json\";\nblah\n";
        jcode.close();
        
        // force the parse to use the path
        wpath = (path / "dump-pretty2.json").native() + "::" + path.native() + ":" + psub.native();
        debug("using WIRECELL_PATH={}", wpath);
        setenv("WIRECELL_PATH", wpath.c_str(), 1);

        auto jroot4 = Persist::load(jpath.native());
        REQUIRE(jroot == jroot4);

    }
    REQUIRE(! boost::filesystem::exists(path));
}
