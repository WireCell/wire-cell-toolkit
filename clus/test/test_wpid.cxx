
#include "WireCellIface/WirePlaneId.h"

#include <iostream>
#include <string>
#include <set>

using namespace WireCell;

int main ()
{
    WireCell::WirePlaneId a5f0pA(kAllLayers, 0, 5);
    WireCell::WirePlaneId a5f1pA(kAllLayers, 1, 5);
    std::cout << "a5f0pA: " << a5f0pA.name() << " ident " << a5f0pA.ident() << std::endl;
    std::cout << "a5f1pA: " << a5f1pA.name() << " ident " << a5f1pA.ident() << std::endl;
    std::cout << " a5f0pA < a5f1pA? " << (a5f0pA < a5f1pA) << std::endl;
    std::cout << " a5f1pA < a5f0pA? " << (a5f1pA < a5f0pA) << std::endl;
    std::cout << " a5f0pA == a5f1pA? " << (a5f0pA == a5f1pA) << std::endl;
    std::cout << " a5f0pA != a5f1pA? " << (a5f0pA != a5f1pA) << std::endl;

    {
        std::set<WireCell::WirePlaneId> wpid_set;
        wpid_set.insert(a5f0pA);
        wpid_set.insert(a5f1pA);
        for (auto wpid : wpid_set) {
            std::cout << "wpid_set [0, 1]: " << wpid.name() << std::endl;
        }
    }
    {
        std::set<WireCell::WirePlaneId> wpid_set;
        wpid_set.insert(a5f1pA);
        wpid_set.insert(a5f0pA);
        for (auto wpid : wpid_set) {
            std::cout << "wpid_set [1, 0]: " << wpid.name() << std::endl;
        }
    }
    {
        std::set<int> int_set;
        int_set.insert(a5f1pA.ident());
        int_set.insert(a5f0pA.ident());
        for (auto wpid_ident : int_set) {
            std::cout << "int_set: " << WirePlaneId(wpid_ident).name() << std::endl;
        }
    }
    return 0;
}