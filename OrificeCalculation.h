#pragma once

class OrificeCalculation
{
  public:
    OrificeCalculation(){}
    ~OrificeCalculation(){}
    
    OrificeCalculation
    (
        double viscosity /* (Pa.s)*/,
        double density /* (kg/m^3)*/,
		double error_tolerance
    );
    
	enum Tapping_Option
	{
		corner_tapping,
		d_and_d_by_2_tapping,
		flange_tapping
	};

    double GetOrificeDiameter
    (     
		double upstream_pressure,	// Pa
		double downstream_pressure, // Pa
		double flow_rate,			// kg/s
		double pipe_in_diameter,	//m
		double isentropic_exponent, //dimensionless
		int tapping_option    
	);
    
    double GetOrificeMassFlowRate
    (
		double upstream_pressure,	// Pa
		double downstream_pressure, // Pa
        double pipe_in_diameter,    // m
		double orifice_diameter	,	//m
		double isentropic_exponent,	// dimensionless
		int tapping_option	    
	);

  private:
  
	double GetOrificeMassFlowRate_NotChoked();

	double GetOrificeMassFlowRate_Choked();

	double GetOrificeDiameter_NotChoked();

	double GetOrificeDiameter_Choked();

	double CalculateDischargeCoefficient(	double beta, 
											double pipe_inner_dia, /* meters*/
											double Reynolds_No_D, 
											Tapping_Option tapping_option);

	double GetReynolds_D_Assumption();

	double CalculateEpsilon(double beta, double pressure_ratio, double isentropic_exponent);

	void GetPressureTappingSpacing(Tapping_Option tapping_option, 
									double pipe_inner_dia, 
									double & L1, 
									double & L2);

    double m_Viscosity;
    double m_Density;
    double m_D;
    double m_d;
	double m_K;
	double m_UpstreamP;
	double m_DownstreamP;
	Tapping_Option m_TappingOption;
	double m_FlowRate;
	double m_ErrorTolerance;
	double m_ReynoldsNo_D;
    //double m_beta;
};