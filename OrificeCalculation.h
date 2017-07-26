#pragma once

class OrificeCalculation
{
  public:
    OrificeCalculation(){}
    ~OrificeCalculation(){}
    
    OrificeCalculation
    (
        double, //viscosity (Pa.s)
        double  //density (kg/m^3)
    );
    
    double diameter
    (
        double pressure_drop, // Pa
        double flow_rate      // kg/s
    );
    
    double flow_rate
    (
        double pressure_drop, // Pa
        double diameter       // m
    );
  
  private:
  
    double viscosity1;
    double density1;
    double D;
    double d;
    double beta;
};