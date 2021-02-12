#ifndef TIMING_HPP
#define TIMING_HPP

#include <chrono>
#include <cstdio>

#define TIME_PT(t)                                       \
    std::chrono::time_point<std::chrono::steady_clock> t \
        = std::chrono::steady_clock::now()

#define TIME_MSG(some_time_1, ...)                     \
    do                                                 \
    {                                                  \
        TIME_PT(some_time_2);                          \
        std::chrono::duration<double> diff             \
            = some_time_2 - some_time_1;               \
        std::fprintf(stderr, "\tTook %.4f sec for ",   \
                             diff.count());            \
        std::fprintf(stderr, __VA_ARGS__);             \
        std::fprintf(stderr, "\n");                    \
    } while (0)

#endif
