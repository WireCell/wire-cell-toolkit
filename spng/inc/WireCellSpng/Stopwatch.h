#include <chrono> // Required for std::chrono::high_resolution_clock, duration, etc.
#include <iostream> // Required for console output (for example usage)

/**
 * @brief A simple stopwatch class using std::chrono::high_resolution_clock.
 *
 * This class provides basic stopwatch functionality, allowing you to start a timer,
 * and then get the elapsed time (lap) in microseconds as a double-precision number.
 */
class Stopwatch {
private:
    // Member variable to store the start time point.
    // high_resolution_clock is generally the best choice for precision.
    std::chrono::time_point<std::chrono::high_resolution_clock> m_startTime;

public:
    /**
     * @brief Constructor for the Stopwatch class.
     *
     * Initializes the stopwatch by calling start() to set an initial time point.
     */
    Stopwatch() {
        start(); // Initialize the start time upon creation
    }

    /**
     * @brief Starts or resets the stopwatch.
     *
     * This method records the current time using high_resolution_clock::now()
     * and stores it as the new starting point for duration calculations.
     * If called again, it effectively resets the timer.
     */
    void start() {
        m_startTime = std::chrono::high_resolution_clock::now();
    }

    /**
     * @brief Restart the stopwatch, return prior duration.
     */
    double restart() {
        const double ret = lap();
        m_startTime = std::chrono::high_resolution_clock::now();
        return ret;
    }

    /**
     * @brief Calculates the elapsed time since the last call to start().
     *
     * @return The duration in microseconds as a double-precision floating-point number.
     * This allows for fractional microseconds if the clock resolution supports it.
     */
    double lap() const {
        // Get the current time point.
        auto currentTime = std::chrono::high_resolution_clock::now();

        // Calculate the duration between the start time and the current time.
        // The duration_cast is used to convert the duration to the desired unit
        // (microseconds) and type (double).
        // std::chrono::duration<double, std::micro> specifies a duration type
        // where the underlying representation is a double and the period is micro (10^-6 seconds).
        std::chrono::duration<double, std::micro> elapsed = currentTime - m_startTime;

        // Return the count of the elapsed duration in microseconds.
        return elapsed.count();
    }
};
