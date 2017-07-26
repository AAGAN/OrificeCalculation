#include "OrificeCalculation.h"
#include <iostream>

int main()
{
    
    OrificeCalculation test_calculation(1e-6,1e-3);
    
    std::cout << test_calculation.diameter(1000.0,1.0) << std::endl;
    std::cout << test_calculation.flow_rate(1000.0,0.01) << std::endl;
    
    return 0;
}