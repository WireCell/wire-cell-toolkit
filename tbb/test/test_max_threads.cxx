#include <tbb/info.h>
#include <tbb/global_control.h>

#include <iostream>

int main()
{
    std::cerr << tbb::info::default_concurrency() << std::endl;
    {
        tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, 2);
        std::cerr << global_limit.active_value(tbb::global_control::max_allowed_parallelism) << std::endl;
        {
            tbb::global_control global_limit2(tbb::global_control::max_allowed_parallelism, 4);
            std::cerr << "\t" <<
                global_limit2.active_value(tbb::global_control::max_allowed_parallelism) << std::endl;
        }
    }

    {
        tbb::global_control* gl = new tbb::global_control(tbb::global_control::max_allowed_parallelism, 16);
        std::cerr << gl->active_value(tbb::global_control::max_allowed_parallelism) << std::endl;
        delete gl;
    }

    return 0;
}
