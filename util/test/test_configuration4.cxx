#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Testing.h"

#include <map>
#include <vector>
#include <iostream>

using namespace WireCell;
int main()
{
    {
        Configuration cfg;
        std::vector<size_t> sizes = {2,4,8};
        assign(cfg, sizes);
        Assert(cfg.size() == 3);
        for (size_t ind=0; ind<sizes.size(); ++ind) {
            auto c = cfg[(int)ind];
            std::cerr << ind << " " << c << "\n";
            Assert(c.asUInt64() == sizes[ind]);
        }
    }
    {
        Configuration cfg;
        std::map<std::string, size_t> named = {{"two",2},{"four",4},{"eight",8}};
        assign(cfg, named);
        Assert(cfg.size() == 3);
        for (const auto& [nam,siz] : named) {
            auto c = cfg[nam];
            std::cerr << nam << " " << c << "\n";
            Assert(c.asUInt64() == siz);
        }
    }
    {
        Configuration cfg;
        std::vector<size_t> sizes = {2,4,8};
        std::map<std::string, std::vector<size_t>> nameds = {{"foo", sizes}, {"bar", {1,2,3}}};
        assign(cfg, nameds);
        Assert(cfg.size() == 2);
        for (const auto& [nam,sizs] : nameds) {
            auto cs = cfg[nam];
            for (size_t ind=0; ind<sizes.size(); ++ind) {
                auto c = cs[(int)ind];
                std::cerr << nam << " " << ind << " " << c << "\n";
                Assert(c.asUInt64() == sizs[ind]);
            }
        }
    }    

    return 0;
}
