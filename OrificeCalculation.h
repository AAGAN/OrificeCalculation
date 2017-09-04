#pragma once

// Error codes
#define ERCODE_ALL_OK 0
#define ERCODE_INVALID_INPUT -1
#define ERCODE_MAX_ITERATIONS_EXCEED -2
#define ERCODE_ERROR_DIVERGING -3
#define ERCODE_PRESSURE_DROP_NEGATIVE -4
#define ERCODE_MACH_NUMBER_NEGATIVE -5
#define ERCODE_REYNOLDS_NUMBER_NEGATIVE -6
#define ERCODE_DISCHARGE_COEFFICIENT_NEGATIVE -7

#define MAX_NO_OF_ITERATIONS 10000
#define MAX_ERROR_CONSECUTIVE_DIVERGE_ALLOWED 20 // we will allow 20 consecutive iterations
												 // where error diverges till we quit
class OrificeCalculation
{
  public:
    OrificeCalculation(){}
    ~OrificeCalculation(){}
    
    OrificeCalculation
    (
        double viscosity	/* (Pa.s)*/,
        double density		/* (kg/m^3)*/,
		double error_tolerance,
		bool compressible,	/* is the fluid compressible?*/
		double gas_const,	/* J/KgK */
		double inlet_temp	/* deg kelvin*/
    );
    
	enum Tapping_Option
	{
		Not_Applicable,
		Corner_Tapping,
		D_And_D_by_2_Tapping,
		Flange_Tapping
	};

    int GetOrificeDiameter
    (     
		double upstream_pressure,	// Pa
		double downstream_pressure, // Pa
		double flow_rate,			// kg/s
		double pipe_in_diameter,	//m
		double isentropic_exponent, //dimensionless
		int tapping_option ,		// enum Tapping_Option values
		double & orifice_dia		// orifice diameter in meters
	);
    
    int GetOrificeMassFlowRate
    (
		double upstream_pressure,	// Pa
		double downstream_pressure, // Pa
        double pipe_in_diameter,    // m
		double orifice_diameter	,	//m
		double isentropic_exponent,	// dimensionless
		int tapping_option,			// enum Tapping_Option values
		double &mass_flow_rate		// 
	);

  private:

	int ValidateInputs();
  
	//int GetOrificeMassFlowRate_InCompressible(double& mass_flow_rate);

	int GetOrificeMassFlowRate_I_5167(double& mass_flow_rate);

	int GetOrificeMassFlowRate_Compressible_Choked(double& mass_flow_rate);

	//int GetOrificeMassFlowRate_Compressible_NotChoked(double& mass_flow_rate);

	//int GetOrificeDiameter_InCompressible(double& mass_flow_rate);

	int GetOrificeDiameter_I_5167(double &diameter);

	int GetOrificeDiameter_Compressible_Choked(double &diameter);

	//int GetOrificeDiameter_Compressible_NotChoked(double& mass_flow_rate);

	int CalculateDischargeCoefficient(		double beta, 
											double pipe_inner_dia, /* meters*/
											double Reynolds_No_D, 
											Tapping_Option tapping_option,
											double & disch_coeff);

	double GetReynolds_D_Assumption();

	double CalculateEpsilon(	double beta, 
								double pressure_ratio, 
								double isentropic_exponent);

	int CalculateMachNumber(	double inlet_pressure, 
								double outlet_pressure,
								double isentropic_coeff, 
								double& Mach_number);

	void GetPressureTappingSpacing(	Tapping_Option tapping_option, 
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
	bool m_Compressible;
	double m_GasConst_R;
	double m_InletTemp;
    //double m_beta;
}; 