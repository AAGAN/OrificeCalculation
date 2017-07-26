#include "OrificeCalculation.h"

OrificeCalculation::OrificeCalculation
    (
        double, //viscosity (Pa.s)
        double  //density (kg/m^3)
    )
    {
        
        
    }
    
double OrificeCalculation::diameter
    (
        double pressure_drop, // Pa
        double flow_rate      // kg/s
    )
    {
        
        return 1.0;
    }

double OrificeCalculation::flow_rate
    (
        double pressure_drop, // Pa
        double diameter       // m
    )
    {
        
        return 1.0;
    }