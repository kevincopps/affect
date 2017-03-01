#ifndef AFFECT_DURATION_H
#define AFFECT_DURATION_H

#include <ctime>
#include <iostream>

#ifdef _DURATION

#define START_DURATION(N)        \
  clock_t start##N, finish##N;  \
  double duration##N;          \
  start##N = clock();
  
#define STOP_DURATION(N,MESSAGE) \
  finish##N = clock();     \
  duration##N = (double)(finish##N - start##N) / CLOCKS_PER_SEC; \
  cout << #MESSAGE << ": " << duration##N << " seconds" << endl;      
  
#else

#define START_DURATION(N)
#define STOP_DURATION(N,STR)
  
#endif

    
#endif // AFFECT_DURATION_H
