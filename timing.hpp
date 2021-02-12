#ifndef TIMING_HPP
#define TIMING_HPP

#include <chrono>
#include <cstdio>

#define TIME_PT(t)                                       \
    std::chrono::time_point<std::chrono::steady_clock> t \
        = std::chrono::steady_clock::now()

#define TIME_MSG(t1, ...)                              \
    do                                                 \
    {                                                  \
        TIME_PT(t2);                                   \
        std::chrono::duration<double> diff = t2 - t1;  \
        std::fprintf(stderr, "\tTook %.4f sec for ",   \
                             diff.count());            \
        std::fprintf(stderr, __VA_ARGS__);             \
        std::fprintf(stderr, "\n");                    \
    } while (0)

#endif
