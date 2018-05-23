#ifndef AFFECT_UTIL_TIMER_H
#define AFFECT_UTIL_TIMER_H

#include <chrono>
#include <iostream>
#include <iomanip>

#ifdef TIMERS_ACTIVE

#define START_TIMER(NAME) \
    auto start##NAME = std::chrono::steady_clock::now();

#define END_TIMER(NAME) \
    auto end##NAME = std::chrono::steady_clock::now(); \
    auto milliseconds##NAME = std::chrono::duration_cast<std::chrono::milliseconds>(end##NAME - start##NAME).count(); \
    std::cout << milliseconds##NAME << "ms " << #NAME << std::endl << std::flush;


#else

#define START_TIMER(NAME)
#define END_TIMER(NAME)

#endif

    //std::chrono::duration<double> elapsed_milliseconds##NAME = end##NAME - start##NAME;
    //std::cerr << std::fixed << std::setprecision(4) << elapsed_seconds##NAME.count() << "s " << #NAME << std::endl << std::flush;

#endif // AFFECT_UTIL_TIMER_H
