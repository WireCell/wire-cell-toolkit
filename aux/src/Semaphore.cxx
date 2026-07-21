#include "WireCellAux/Semaphore.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Semaphore.h"

WIRECELL_FACTORY(Semaphore, 
                 WireCell::Aux::Semaphore,
                 WireCell::ISemaphore,
                 WireCell::IConfigurable)

using namespace WireCell;

Aux::Semaphore::Semaphore()
    : m_sem(1)
{
}
Aux::Semaphore::~Semaphore()
{
}

WireCell::Configuration Aux::Semaphore::default_configuration() const
{
    Configuration cfg;

    // The maximum allowed number concurrent calls to forward().  A
    // value of unity means all calls will be serialized.  When made
    // smaller than the number of threads, the difference gives the
    // number of threads that may block on the semaphore.
    cfg["concurrency"] = 1;

    return cfg;
}

void Aux::Semaphore::configure(const WireCell::Configuration& cfg)
{
    auto count = get<int>(cfg, "concurrency", 1);
    if (count < 1 ) {
        count = 1;
    }
    m_sem.set_count(count);
}

void Aux::Semaphore::acquire() const
{
    m_sem.acquire();
}

void Aux::Semaphore::release() const
{
    m_sem.release();
}

int Aux::Semaphore::get_count() const
{
    m_sem.get_count();
}
