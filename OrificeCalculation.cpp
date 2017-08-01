#include "Math.h"
#include "OrificeCalculation.h"

#define CHOCKING_PRESSURE_RATIO 0.75 // required?
#define ATMOSPHERIC_PRESSURE_Pa 101325 // required?
#define _PI_ 3.14159265359
#define _METER_TO_MM_ 1000
#define _e_ 2.718281828

OrificeCalculation::OrificeCalculation
    (
        double viscosity /* (Pa.s) */,
        double  density /*(kg/m^3) */,
		double error_tolerance
    )
    {
		m_Viscosity = viscosity;
		m_Density = density;
		m_ErrorTolerance = error_tolerance;
    }
    
	double OrificeCalculation::GetOrificeDiameter (double upstream_pressure, double downstream_pressure,
												  double flow_rate, double pipe_in_diameter,
												  double isentropic_exponent, int tapping_option)
    {
		double result = 0;
		m_UpstreamP = upstream_pressure;
		m_DownstreamP = downstream_pressure;
		m_D = pipe_in_diameter;
		m_FlowRate = flow_rate;
		m_K = isentropic_exponent;
		m_TappingOption = (Tapping_Option) tapping_option;
		double pressure_ratio = downstream_pressure / upstream_pressure;
		if (pressure_ratio > CHOCKING_PRESSURE_RATIO)
			result = GetOrificeDiameter_NotChoked();
		else
			result = GetOrificeDiameter_Choked();
        return result;
    }

	double OrificeCalculation::GetOrificeDiameter_NotChoked()
	{
		double result = 0;
		double X2_final = 0;
		// Calculate Reynolds No for Diameter
		double Re_D = 4 * m_FlowRate / (_PI_ * m_Viscosity * m_D);

		// Calculate Invariant
		double pressure_diff = m_UpstreamP - m_DownstreamP;
		double pressure_ratio = m_DownstreamP / m_UpstreamP;
		double A2 = m_Viscosity * Re_D / (m_D * sqrt(2 * pressure_diff * m_Density));

		bool error_larger_than_tolerance = true;
		double beta_1 = 0.2;		// first assumption
		double C1 = CalculateDischargeCoefficient(beta_1, m_D, Re_D, m_TappingOption);
		double epsilon1 = CalculateEpsilon(beta_1, pressure_ratio, m_K);
		double X2_1 = pow(beta_1, 2) / sqrt(1 - pow(beta_1, 4));
		double delta1 = fabs( (A2 - X2_1* C1* epsilon1) / A2);

		if (delta1 < m_ErrorTolerance)
		{
			error_larger_than_tolerance = false;
			X2_final = X2_1;
		}
		double beta_2 = 0.1;	// second assumption
		double C2 = CalculateDischargeCoefficient(beta_2, m_D, Re_D, m_TappingOption);
		double epsilon2 = CalculateEpsilon(beta_2, pressure_ratio, m_K);
		double X2_2 = pow(beta_2, 2) / sqrt(1 - pow(beta_2, 4));
		double delta2 = fabs((A2 - X2_2* C2* epsilon2) / A2);

		if (error_larger_than_tolerance && delta2 < m_ErrorTolerance)
		{
			error_larger_than_tolerance = false;
			X2_final = X2_2;
		}

		double X2_1_before = X2_1;
		double X2_2_before = X2_2;
		double er_1_before = delta1;
		double er_2_before = delta2;
		double beta_temp = 0;
		while (error_larger_than_tolerance)
		{	// some check to see if fails to converge?
			double X2_temp = X2_1_before - er_1_before * (X2_1_before - X2_2_before) / (er_1_before - er_2_before);
			double var_temp = X2_temp* X2_temp / (1 + X2_temp* X2_temp);
			double beta_temp = pow(  X2_temp* X2_temp / (1 + X2_temp* X2_temp) , 0.25);
			double C_temp = CalculateDischargeCoefficient(beta_temp, m_D, Re_D, m_TappingOption);
			double epsilon_temp = CalculateEpsilon(beta_temp, pressure_ratio, m_K);

			double error = fabs((A2 - X2_temp * C_temp* epsilon_temp) / A2) ;
			if (error < m_ErrorTolerance)
			{
				error_larger_than_tolerance = false;
				X2_final = X2_temp;
				break;
			}
			else
			{
				X2_1_before = X2_2_before;
				X2_2_before = X2_temp;
				er_1_before = er_2_before;
				er_2_before = error;
				continue;
			}
		}

		double beta_final = pow((X2_final* X2_final / (1 + X2_final* X2_final)), 0.25);
		result = beta_final * m_D;

		return result;
	}

	double OrificeCalculation::GetOrificeDiameter_Choked()
	{
		double result = 0;
		return result;
	}

	double OrificeCalculation::GetOrificeMassFlowRate ( double upstream_pressure, double downstream_pressure,
														double pipe_in_diameter,  double orifice_diameter,
														double isentropic_exponent, int tapping_option)
    {
		double result = 0.0;
		m_UpstreamP = upstream_pressure;
		m_DownstreamP = downstream_pressure;
		m_D = pipe_in_diameter;	
		m_d = orifice_diameter;	
		m_K = isentropic_exponent;
		m_TappingOption = (Tapping_Option) tapping_option;
		double pressure_ratio = downstream_pressure / upstream_pressure;
		if (pressure_ratio > CHOCKING_PRESSURE_RATIO)
			result = GetOrificeMassFlowRate_NotChoked();
		else
			result = GetOrificeMassFlowRate_Choked();
        return result;
    }

	double OrificeCalculation::GetOrificeMassFlowRate_NotChoked()
	{
		double result = 0;

		// Calculate Beta
		double beta = m_d / m_D;

		// Calculate expansibility factor, epsilon
		double pressure_ratio = m_DownstreamP / m_UpstreamP;
		double pressure_drop = m_UpstreamP - m_DownstreamP;
		double epsilon = CalculateEpsilon(beta, pressure_ratio, m_K); 

		// Calculate the invariant
		double A1_numerator = epsilon * pow(m_d, 2) * sqrt(2 * pressure_drop * m_Density);
		double A1_denominator = m_Viscosity * m_D * sqrt(1 - pow(beta, 4));
		double A1 = A1_numerator / A1_denominator;

		// Get two assumptions of Reynolds Number D
		double X1 = GetReynolds_D_Assumption() * 1.1;
		double X2 = GetReynolds_D_Assumption() * 1.4;

		// Calculate discharge coefficients C1,C2
		double C1 = CalculateDischargeCoefficient(beta, m_D, X1, m_TappingOption);
		double C2 = CalculateDischargeCoefficient(beta, m_D, X2, m_TappingOption);

		//calculate two errors
		bool error_larger_than_tolerance = true;
		double X= 0;
		double delta1 = fabs((A1 - X1 / C1) / A1);
		if (delta1 < m_ErrorTolerance)
		{ 
			error_larger_than_tolerance = false;
			X = X1;
		}
		double delta2 = fabs((A1 - X2 / C2) / A1);
		if (error_larger_than_tolerance && delta2 < m_ErrorTolerance)
		{ 
			error_larger_than_tolerance = false;
			X = X2;
		}
		double X_1_before = X1;
		double X_2_before = X2;
		double er_1_before = delta1;
		double er_2_before = delta2;
		while (error_larger_than_tolerance)
		{	// some check to see if fails to converge?
			double X_temp = X_1_before - er_1_before * (X_1_before - X_2_before) / (er_1_before - er_2_before);
			double C_temp = CalculateDischargeCoefficient(beta, m_D, X_temp, m_TappingOption);
			double error = fabs((A1 - X_temp / C_temp) / A1);
			if (error < m_ErrorTolerance)
			{
				error_larger_than_tolerance = false;
				X = X_temp;
				break;
			}
			else
			{
				X_1_before = X_2_before;
				X_2_before = X_temp;
				er_1_before = er_2_before;
				er_2_before = error;
				continue;
			}
		}

		result = _PI_ / 4 * m_Viscosity * m_D * X;
		return result;

	}

	double OrificeCalculation::GetOrificeMassFlowRate_Choked()
	{
		double result = 0;
		return result;
	}

	double OrificeCalculation::CalculateEpsilon(double beta, double pressure_ratio, double isentropic_exponent)
	{
		// ISO 5167-2:2003(E) , Section 5.2.2.2
		double epsilon = 1 - (0.351 + 0.256 * pow(beta, 4) + 0.93 * pow(beta, 8)) * (1 - pow(pressure_ratio, (1 / isentropic_exponent)));
		return epsilon;
	}

	void OrificeCalculation::GetPressureTappingSpacing(Tapping_Option tapping_option, double pipe_inner_dia, double & L1, double & L2)
	{
		switch (tapping_option) // make this a function
		{
			case (Tapping_Option::corner_tapping):
			{
				L1 = 0;
				L2 = 0;
			}
			break;
			case (Tapping_Option::d_and_d_by_2_tapping):
			{
				L1 = 1;
				L2 = 0.47;
			}
			break;
			case (Tapping_Option::flange_tapping):
			{
				L1 = 25.4 / (pipe_inner_dia* _METER_TO_MM_);
				L2 = L1;
			}
			break;
		}
	}

	double OrificeCalculation::CalculateDischargeCoefficient( double beta, double pipe_inner_dia, 
															  double Reynolds_No_D, Tapping_Option tapping_option)
	{
		// Formula from ISO 5167-2:2003(E) section 5.3.2.1
		//double beta = m_d / m_D;
		double A = pow((19000 * beta / Reynolds_No_D), 0.8);
		double L1, L2;
		GetPressureTappingSpacing(tapping_option, pipe_inner_dia, L1, L2);
		double M2 = 2 * L2 / (1 - beta);
		double C = 0.5961 +
			0.0261 * pow(beta, 2) - 0.216 * pow(beta, 8) +
			0.000521 * pow((10e6 * beta / Reynolds_No_D), 0.7) +
			(0.0188 + 0.0063 * A) * pow(beta, 3.5) * pow((10e6 / Reynolds_No_D), 0.3) +
			(0.043 + 0.08 * pow(_e_, (-10 * L1)) - 0.123 * pow(_e_, -7 * L1)) * (1 - 0.11 * A) * (pow(beta, 4) / (1 - pow(beta, 4))) -
			0.031 * (M2 - 0.8 *  pow(M2, 1.1) * pow(beta, 1.3));

		double C_adjustment = 0;
		if (m_D * _METER_TO_MM_ < 71.12 )
		{
			C_adjustment = 0.011 * (0.75 - beta) * (2.8 - m_D * _METER_TO_MM_ / 25.4);
		}
		C += C_adjustment;

		return C;

	}

	double OrificeCalculation::GetReynolds_D_Assumption()
	{
		double result = 0;
		double beta = m_d / m_D;
		switch (m_TappingOption)
		{
			case Tapping_Option::corner_tapping:
			case Tapping_Option::d_and_d_by_2_tapping:
			{
				if (beta >= 0.1 && beta <= 0.56)
					result = 5000;
				else if (beta > 0.56)
					result = 16000;
				else
					result = 0; // this should be error?
			}
			break;
			case Tapping_Option::flange_tapping:
			{
				result = 170 * pow(beta, 2) * m_D * _METER_TO_MM_;
				if (result < 5000)
					result = 5000;
			}
			break;
		}
		return result;
	}
