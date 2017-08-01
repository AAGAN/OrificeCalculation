#include "OrificeCalculation.h"
#include <iostream>

int main()
{
    
    OrificeCalculation test_calculation(4.0416e-6,4.0175, 1e-4);
    
    std::cout << test_calculation.GetOrificeMassFlowRate( 600000, 500000, 0.1, 0.025, 1.3, 0) << std::endl;
    //std::cout << test_calculation.GetOrificeDiameter( 600000, 500000, 0.248, 0.1, 1.3, 0) << std::endl;
    
    return 0;
}