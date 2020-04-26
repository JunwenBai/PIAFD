//
//  ChronoP.hpp
//  MDD
//
//  Created by guillaume perez on 04/11/2015.
//  Copyright © 2015 MemoCop. All rights reserved.
//

#ifndef ChronoP_hpp
#define ChronoP_hpp

#include <iostream>
#include <chrono>
#include <ctime>

class ChronoP {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    
public:
    void Start();
    void Stop();
    void Restart();
    
    int64_t ellapsed_second();
    int64_t ellapsed_m_second();
    int64_t ellapsed_u_second();
    
    int64_t compute_granularity();
};

#endif /* ChronoP_hpp */
