/** Test the WireCellUtil/Testing helpers.
 */

#include "WireCellUtil/Testing.h"
#include "WireCellUtil/ExecMon.h"

#include <unistd.h>  // fixme: maybe a more modern sleep() is preferable?
#include <iostream>

using namespace WireCell;
using namespace std;

int fact(int n) {
    if (n <= 1) {
        return 1;
    }
    return n*fact(n-1);
}

void spin_piggy()
{
    int big = 100;    
    double dummy=0;
    while (big) {
        --big;
        dummy += fact(big);
    }
}

void spin_piggy_em()
{
    ExecMon em("");
    int big = 100;    
    double dummy=0;
    while (big) {
        --big;
        dummy += fact(big);
        em("");
    }
    cout << em("last") << endl;
}


void test_speed()
{
    
}

int main(int argc, char** argv)
{
    ExecMon em("test_testing");
    cout << em("sleeping") << endl;
    ;
    sleep(1);
    cout << em("awake") << endl;

    spin_piggy();
    cout << em("no em") << endl;
    spin_piggy();
    cout << em("yes em") << endl;

    cout << em.summary() << endl;
    return 0;
}
