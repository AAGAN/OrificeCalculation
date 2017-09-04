#include "Math.h"
#include "OrificeCalculation.h"

#define _I_5167_PRESSURE_RATIO_LIMIT_ 0.75 // required?
#define _CHOKING_PRESSURE_RATIO_UPPER_LIMIT_ 0.587
#define _CHOKING_PRESSURE_RATIO_LOWER_LIMIT_ 0.487

#define ATMOSPHERIC_PRESSURE_Pa 101325 // required?
#define _PI_ 3.14159265359 
#define _METER_TO_MM_ 1000
#define _e_ 2.718281828
#define _ORIFICE_FLOW_COEFF_ 0.6

OrificeCalculation::OrificeCalculation
    (
        double viscosity /* (Pa.s) */,
        double  density /*(kg/m^3) */,
		double error_tolerance,
		bool compressible,
		double gas_const,
		double inlet_temp
    )
    {
		m_Viscosity = viscosity;
		m_Density = density;
		m_ErrorTolerance = error_tolerance;
		m_Compressible = compressible;
		m_GasConst_R = gas_const;
		m_InletTemp = inlet_temp;
    }

	int OrificeCalculation::ValidateInputs()
	{
		if (m_Viscosity < 0 || m_Density < 0 || m_Density < 0 || m_GasConst_R < 0 || m_InletTemp < 0
			|| m_UpstreamP < 0 || m_DownstreamP < 0 || m_D < 0
			|| (m_d < 0 && m_FlowRate < 0) || m_K < 0
			|| m_TappingOption < 0)
			return ERCODE_INVALID_INPUT;
		else if (m_UpstreamP < m_DownstreamP)
			return ERCODE_PRESSURE_DROP_NEGATIVE;
		else
			return ERCODE_ALL_OK;
	}
    
	int OrificeCalculation::GetOrificeDiameter (double upstream_pressure, double downstream_pressure,
												  double flow_rate, double pipe_in_diameter,
												  double isentropic_exponent, int tapping_option,
													double &orifice_dia)
    {
		int err_code = 0;
		double result = 0;
		m_UpstreamP = upstream_pressure;
		m_DownstreamP = downstream_pressure;
		m_D = pipe_in_diameter;
		m_FlowRate = flow_rate;
		m_K = isentropic_exponent;
		m_TappingOption = (Tapping_Option) tapping_option;

		err_code = ValidateInputs();
		if (err_code < 0)
			return err_code;
		if (m_Compressible == false)
		{
			err_code = GetOrificeDiameter_I_5167(orifice_dia);
		}
		else
		{
			double pressure_ratio = downstream_pressure / upstream_pressure;
			if (pressure_ratio > _I_5167_PRESSURE_RATIO_LIMIT_)
				err_code = GetOrificeDiameter_I_5167(orifice_dia);
			else
				err_code = GetOrificeDiameter_Compressible_Choked(orifice_dia);
		}

		return err_code;
    }
/*
	int OrificeCalculation::GetOrificeDiameter_InCompressible(double& orifice_dia)
	{
		int er_code = 0;
		double orifice_area = m_FlowRate / (_ORIFICE_FLOW_COEFF_ * sqrt(2 * m_Density * (m_UpstreamP - m_DownstreamP)));
		orifice_dia = sqrt(orifice_area * 4 / _PI_);
		return er_code;
	}
*/

	int OrificeCalculation::GetOrificeDiameter_I_5167(double& diameter)
	{
		int err_code = 0;
		double X2_final = 0;
		// Calculate Reynolds No for Diameter
		double Re_D = 4 * m_FlowRate / (_PI_ * m_Viscosity * m_D);

		// Calculate Invariant
		double pressure_diff = m_UpstreamP - m_DownstreamP;
		double pressure_ratio = m_DownstreamP / m_UpstreamP;
		double A2 = m_Viscosity * Re_D / (m_D * sqrt(2 * pressure_diff * m_Density));

		bool error_larger_than_tolerance = true;
		double beta_1 = 0.2;		// first assumption
		double C1 = 0;
		err_code = CalculateDischargeCoefficient(beta_1, m_D, Re_D, m_TappingOption, C1);
		if (err_code < 0)
			return err_code;
		double epsilon1 = CalculateEpsilon(beta_1, pressure_ratio, m_K);
		double X2_1 = pow(beta_1, 2) / sqrt(1 - pow(beta_1, 4));
		double delta1 = fabs( (A2 - X2_1* C1* epsilon1) / A2);

		if (delta1 < m_ErrorTolerance)
		{
			error_larger_than_tolerance = false;
			X2_final = X2_1;
		}
		double beta_2 = 0.1;	// second assumption
		double C2 = 0;
		err_code = CalculateDischargeCoefficient(beta_2, m_D, Re_D, m_TappingOption,C2);
		if (err_code < 0)
			return err_code;
		double epsilon2 = CalculateEpsilon(beta_2, pressure_ratio, m_K);
		double X2_2 = pow(beta_2, 2) / sqrt(1 - pow(beta_2, 4));
		double delta2 = fabs((A2 - X2_2* C2* epsilon2) / A2);

		if (error_larger_than_tolerance && delta2 < m_ErrorTolerance)
		{
			error_larger_than_tolerance = false;
			X2_final = X2_2;
		}

		// Start iterating till error stabilizes
		double X2_1_before = X2_1;		// X2 before 1 generation
		double X2_2_before = X2_2;		// X2 before 2 generation
		double er_1_before = delta1;	// error before 1 generation
		double er_2_before = delta2;	// error before 2 generation
		double beta_temp = 0;
		int no_of_iterations = 0;
		int err_converge_tracker = 0;
		while (error_larger_than_tolerance)
		{	// some check to see if fails to converge?
			++no_of_iterations;
			if (no_of_iterations > MAX_NO_OF_ITERATIONS)
				return ERCODE_MAX_ITERATIONS_EXCEED;
			if (er_2_before > er_1_before)
				err_converge_tracker += 1;
			else
				err_converge_tracker = 0;
			if (err_converge_tracker > MAX_ERROR_CONSECUTIVE_DIVERGE_ALLOWED)	// not converged for last six iterations
				return ERCODE_ERROR_DIVERGING;
			double X2_temp = X2_1_before - er_1_before * (X2_1_before - X2_2_before) / (er_1_before - er_2_before);
			double var_temp = X2_temp* X2_temp / (1 + X2_temp* X2_temp);
			double beta_temp = pow(  X2_temp* X2_temp / (1 + X2_temp* X2_temp) , 0.25);
			double C_temp = 0;
			err_code = CalculateDischargeCoefficient(beta_temp, m_D, Re_D, m_TappingOption, C_temp);
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
		diameter = beta_final * m_D;

		return err_code;
	}

	int OrificeCalculation::GetOrificeDiameter_Compressible_Choked(double& diameter)
	{
		// The formula here is from
		// http://engineeringstudymaterial.net/ebook/fluid-mechanics-by-yunus-cengel-john/
		// chapter 12		
		int err_code = 0;

		double temp_calc_1 = sqrt(m_K / (m_GasConst_R* m_InletTemp));
		double temp_calc_2 = pow(
			(2 / (m_K + 1)),
			(m_K + 1) / (2 * (m_K - 1))
		);
		double Area_throat = m_FlowRate/  (m_UpstreamP * temp_calc_1 * temp_calc_2);
		diameter = sqrt(4 * Area_throat / _PI_);
		return err_code;
	}
	/*
	int OrificeCalculation::GetOrificeDiameter_Compressible_NotChoked(double& diameter)
	{
		int err_code = 0;
		double temp_calc_1 = m_UpstreamP / sqrt(m_InletTemp) * sqrt(m_K /  m_GasConst_R);
		double mach_number = 0;
		err_code = CalculateMachNumber(m_UpstreamP, m_DownstreamP, m_K, mach_number);
		double temp_calc_2 = mach_number * (1 + (m_K - 1) / 2 * pow(mach_number, 2) );
		double temp_calc_3 = -1 * (m_K + 1) / (2 * (m_K - 1));
		double orifice_area = m_FlowRate / (temp_calc_1 * pow(temp_calc_2, temp_calc_3));
		diameter = sqrt( 4 * orifice_area/ _PI_ );
		return err_code;
	}
	*/
	int OrificeCalculation::GetOrificeMassFlowRate ( double upstream_pressure, double downstream_pressure,
														double pipe_in_diameter,  double orifice_diameter,
														double isentropic_exponent, int tapping_option,
														double &flow_rate)
    {
		int err_code = 0;
		m_UpstreamP = upstream_pressure;
		m_DownstreamP = downstream_pressure;
		m_D = pipe_in_diameter;	
		m_d = orifice_diameter;	
		m_K = isentropic_exponent;
		m_TappingOption = (Tapping_Option) tapping_option;

		// is the fluid compressible or not?
		if (m_Compressible == false) // incompressible
		{
			err_code = GetOrificeMassFlowRate_I_5167(flow_rate);
		}
		else
		{
			double pressure_ratio = downstream_pressure / upstream_pressure;
			if (pressure_ratio > _I_5167_PRESSURE_RATIO_LIMIT_)
				err_code = GetOrificeMassFlowRate_I_5167(flow_rate);
			else
				err_code = GetOrificeMassFlowRate_Compressible_Choked(flow_rate);
		}
		return err_code;
    }
/*
	int OrificeCalculation::GetOrificeMassFlowRate_InCompressible(double& flow_rate)
	{
		int er_code = 0;
		flow_rate = _ORIFICE_FLOW_COEFF_ *
					(_PI_ *  pow(m_d, 2) / 4) *
					sqrt(2 * m_Density * (m_UpstreamP - m_DownstreamP));
		return er_code;
	}
*/
	int OrificeCalculation::GetOrificeMassFlowRate_I_5167(double & flow_rate)
	{
		int err_code = 0;

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
		double C1, C2 = 0;
		err_code = CalculateDischargeCoefficient(beta, m_D, X1, m_TappingOption, C1);
		err_code = CalculateDischargeCoefficient(beta, m_D, X2, m_TappingOption, C2);
		if (C1 < 0 || C2 < 0)
			return ERCODE_DISCHARGE_COEFFICIENT_NEGATIVE;

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
		int no_of_iterations = 0;
		int err_is_converging = 1;
		while (error_larger_than_tolerance)
		{	// some check to see if fails to converge?
			++no_of_iterations;
			if (no_of_iterations > MAX_NO_OF_ITERATIONS)
				return ERCODE_MAX_ITERATIONS_EXCEED;
			if (er_2_before > er_1_before)
				err_is_converging *= 2;
			else
				err_is_converging = 1;
			if (err_is_converging > 64)	// not converged for last six iterations
				return ERCODE_ERROR_DIVERGING;
			double X_temp = X_1_before - er_1_before * (X_1_before - X_2_before) / (er_1_before - er_2_before);
			double C_temp = 0;
			err_code = CalculateDischargeCoefficient(beta, m_D, X_temp, m_TappingOption, C_temp);
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

		flow_rate = _PI_ / 4 * m_Viscosity * m_D * X;
		return err_code;

	}

	int OrificeCalculation::GetOrificeMassFlowRate_Compressible_Choked(double &flow_rate)
	{
		// The formula here is from
		// http://engineeringstudymaterial.net/ebook/fluid-mechanics-by-yunus-cengel-john/
		// chapter 12
		int err_code = 0;
		double temp_calc_1 = sqrt(m_K / (m_GasConst_R* m_InletTemp));
		double temp_calc_2 = pow(
			(2 / (m_K + 1)),
			(m_K + 1) / (2 * (m_K - 1))
		);
		double Area_throat = _PI_ * pow(m_d, 2) / 4;
		flow_rate = Area_throat * m_UpstreamP * temp_calc_1 * temp_calc_2;
		return err_code;
	}
/*
	int OrificeCalculation::GetOrificeMassFlowRate_Compressible_NotChoked(double &flow_rate)
	{
		int err_code = 0;
		double temp_calc_1 = m_UpstreamP / sqrt(m_InletTemp) * sqrt(m_K /  m_GasConst_R);
		double mach_number = 0;
		err_code = CalculateMachNumber(m_UpstreamP, m_DownstreamP, m_K, mach_number);
		double temp_calc_2 = mach_number * (1 + (m_K - 1) / 2 * pow(mach_number, 2) );
		double temp_calc_3 = -1 * (m_K + 1) / (2 * (m_K - 1));
		double orifice_area = _PI_ * pow(m_d,2) / 4;
		flow_rate = orifice_area * temp_calc_1 * pow( temp_calc_2, temp_calc_3);
		return err_code;
	}
*/
	double OrificeCalculation::CalculateEpsilon(double beta, double pressure_ratio, double isentropic_exponent)
	{
		// ISO 5167-2:2003(E) , Section 5.2.2.2
		double epsilon = 0;
		if (!m_Compressible)
			epsilon = 1;
		else
			epsilon = 1 - (0.351 + 0.256 * pow(beta, 4) + 0.93 * pow(beta, 8)) * (1 - pow(pressure_ratio, (1 / isentropic_exponent)));
		return epsilon;
	}

	void OrificeCalculation::GetPressureTappingSpacing(Tapping_Option tapping_option, double pipe_inner_dia, double & L1, double & L2)
	{
		switch (tapping_option) // make this a function
		{
			case (Tapping_Option::Corner_Tapping):
			{
				L1 = 0;
				L2 = 0;
			}
			break;
			case (Tapping_Option::D_And_D_by_2_Tapping):
			{
				L1 = 1;
				L2 = 0.47;
			}
			break;
			case (Tapping_Option::Flange_Tapping):
			{
				L1 = 25.4 / (pipe_inner_dia* _METER_TO_MM_);
				L2 = L1;
			}
			break;
		}
	}

	int OrificeCalculation::CalculateDischargeCoefficient( double beta, double pipe_inner_dia, 
															  double Reynolds_No_D, Tapping_Option tapping_option, 
															double & disch_coeff)
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

		disch_coeff = C;
		if (C < 0)
			return ERCODE_DISCHARGE_COEFFICIENT_NEGATIVE;
		else
			return ERCODE_ALL_OK;

	}

	int OrificeCalculation::CalculateMachNumber(double inlet_pr, double outlet_pr,
		double isentropic_coeff, double& Mach_number)
	{
		int er_code = 0;
		double calc_step_1 = pow((inlet_pr / outlet_pr), (isentropic_coeff - 1) / isentropic_coeff) - 1;
		Mach_number = pow(calc_step_1 * 2 / (isentropic_coeff - 1), 0.5);
		if (Mach_number < 0)
			er_code = ERCODE_MACH_NUMBER_NEGATIVE;
		return er_code;
	}

	double OrificeCalculation::GetReynolds_D_Assumption()
	{
		double result = 0;
		double beta = m_d / m_D;
		switch (m_TappingOption)
		{
			case Tapping_Option::Corner_Tapping:
			case Tapping_Option::D_And_D_by_2_Tapping:
			{
				if (beta >= 0.1 && beta <= 0.56)
					result = 5000;
				else if (beta > 0.56)
					result = 16000;
				else
					result = 0; // this should be error?
			}
			break;
			case Tapping_Option::Flange_Tapping:
			{
				result = 170 * pow(beta, 2) * m_D * _METER_TO_MM_;
				if (result < 5000)
					result = 5000;
			}
			break;
		}
		return result;
	}
