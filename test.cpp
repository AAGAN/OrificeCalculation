#include "OrificeCalculation.h"
#include <iostream>

int main()
{
	// Test case 1: Orifice flow rate
	if (false)
	{
		// Input values
		double viscosity = 4.0416e-6;
		double density = 4.0175;
		double error_tolerance = 1e-4;
		bool compressible = true;
		double gas_constant = 252.502;
		double inlet_temp = 300;

		double upstream_pressure = 600000;
		double downstream_pressure = 500000;
		double pipe_inner_dia = 0.1;
		double orifice_dia = 0.025;
		double isentropic_exponent = 1.3;
		OrificeCalculation::Tapping_Option tapping_option = OrificeCalculation::Tapping_Option::Corner_Tapping;

		// Output values
		int er_code = 0;
		double flow_rate = 0;

		OrificeCalculation test_calculation(viscosity, density,error_tolerance, 
												compressible, gas_constant, inlet_temp);
		er_code = test_calculation.GetOrificeMassFlowRate(upstream_pressure, downstream_pressure, 
														  pipe_inner_dia, orifice_dia, isentropic_exponent,
														  tapping_option, 
														  flow_rate);

		_ASSERT(er_code == 0 && fabs(flow_rate - 0.2520) < 1e-4); // flow_rate should be 0.2520
	}

	// Test case 2 : Orifice diameter
	if (false)
	{
		// Input values
		double viscosity = 4.0416e-6;
		double density = 4.0175;
		double error_tolerance = 1e-4;
		bool compressible = true;
		double gas_constant = 252.502;
		double inlet_temp = 300;

		double upstream_pressure = 600000;
		double downstream_pressure = 500000;
		double pipe_inner_dia = 0.1;
		double orifice_dia = 0.025;
		double isentropic_exponent = 1.3;
		OrificeCalculation::Tapping_Option tapping_option = OrificeCalculation::Tapping_Option::Corner_Tapping;

		// Output values
		int er_code = 0;
		double diameter = 0;

		OrificeCalculation test_calculation(viscosity, density,error_tolerance, compressible, 
																		gas_constant, inlet_temp);	
		er_code = test_calculation.GetOrificeDiameter( 600000, 500000, 0.248, 0.1, 1.3, 0, diameter);
		_ASSERT(er_code == 0 &&  fabs(diameter - 0.02480) < 1e-4); // diameter shold be 0.02480	
	}

	// Test case 3 : Orifice flowrate, incompressible
	if (true)
	{
		// Input values
		double viscosity = 1.006e-6;
		double density = 1000;
		double error_tolerance = 1e-4;
		bool compressible = false;
		double gas_constant = 252.502;
		double inlet_temp = 300;

		double upstream_pressure = 600000;
		double downstream_pressure = 500000;
		double pipe_inner_dia = 0.1;
		double orifice_dia = 0.025;
		double isentropic_exponent = 1.3;
		OrificeCalculation::Tapping_Option tapping_option = OrificeCalculation::Tapping_Option::Corner_Tapping;

		// Output values
		int er_code = 0;
		double flow_rate = 0;

		OrificeCalculation test_calculation(viscosity, density,error_tolerance, compressible, 
																		gas_constant, inlet_temp);
		er_code = test_calculation.GetOrificeMassFlowRate(upstream_pressure, downstream_pressure, 
														  pipe_inner_dia, orifice_dia, isentropic_exponent,
														  tapping_option, 
														  flow_rate);
		_ASSERT(er_code == 0 && fabs(flow_rate - 4.17) < 1e-2); // flow_rate should be 4.17
	}

	//Test case 4: Orifice diameter, incompressible
	if (true)
	{
		// Input values
		double viscosity = 1.006e-6;
		double density = 1000;
		double error_tolerance = 1e-4;
		bool compressible = false;
		double gas_constant = 252.502;
		double inlet_temp = 300;

		double flow_rate = 4.17;
		double upstream_pressure = 600000;
		double downstream_pressure = 500000;
		double pipe_inner_dia = 0.1;
		double isentropic_exponent = 1.3;
		OrificeCalculation::Tapping_Option tapping_option = OrificeCalculation::Tapping_Option::Corner_Tapping;

		// Output values
		int er_code = 0;
		double orifice_dia = 0;

		OrificeCalculation test_calculation(viscosity, density,error_tolerance, compressible, 
																		gas_constant, inlet_temp);
		er_code = test_calculation.GetOrificeDiameter(upstream_pressure, downstream_pressure, 
														  flow_rate, pipe_inner_dia, isentropic_exponent,
														  tapping_option, 
														  flow_rate);
		_ASSERT(er_code == 0 && fabs(orifice_dia - 0.025) < 1e-2); // flow_rate should be 0.025
	}

	// Test case 4 : Orifice flowrate, tapping option - flange
	if(false)
	{
		// Input values
		double viscosity = 4.0416e-6;
		double density = 4.0175;
		double error_tolerance = 1e-4;
		bool compressible = true;
		double gas_constant = 252.502;
		double inlet_temp = 300;

		double upstream_pressure = 600000;
		double downstream_pressure = 500000;
		double pipe_inner_dia = 0.1;
		double orifice_dia = 0.025;
		double isentropic_exponent = 1.3;
		OrificeCalculation::Tapping_Option tapping_option = OrificeCalculation::Tapping_Option::Flange_Tapping;

		// Output values
		int er_code = 0;
		double flow_rate = 0;

		OrificeCalculation test_calculation(viscosity, density,error_tolerance, compressible, 
																		gas_constant, inlet_temp);
		er_code = test_calculation.GetOrificeMassFlowRate(upstream_pressure, downstream_pressure, 
														  pipe_inner_dia, orifice_dia, isentropic_exponent,
														  tapping_option, 
														  flow_rate);

		//_ASSERT(er_code == 0 && fabs(flow_rate - 4.17) < 1e-2); // flow_rate should be 4.17
	}

	// Test case 5 : Choked flow, compressible, calculate mass flow rate
	if(false)
	{
		// Input values
		double viscosity = 4.0416e-6;
		double density = 4.0175;
		double error_tolerance = 1e-4;
		bool compressible = true;
		double gas_constant = 287.1;
		double inlet_temp = 350;

		double upstream_pressure = 1000000;
		double downstream_pressure = 500000;
		double pipe_inner_dia = 0.03568248232;
		double orifice_dia = 0.0248;
		double isentropic_exponent = 1.4;
		OrificeCalculation::Tapping_Option tapping_option = OrificeCalculation::Tapping_Option::Flange_Tapping;

		// Output values
		int er_code = 0;
		double flow_rate = 0;

		OrificeCalculation test_calculation(viscosity, density,error_tolerance, compressible, 
																		gas_constant, inlet_temp);
		er_code = test_calculation.GetOrificeMassFlowRate(upstream_pressure, downstream_pressure, 
														  pipe_inner_dia, orifice_dia, isentropic_exponent,
														  tapping_option, 
														  flow_rate);

		//_ASSERT(er_code == 0 && fabs(flow_rate - 4.17) < 1e-2); // flow_rate should be 4.17
	}

	// Test case 5 : Choked flow, compressible, calculate orifice dia
	if(true)
	{
		// Input values
		double viscosity = 4.0416e-6;
		double density = 4.0175;
		double error_tolerance = 1e-4;
		bool compressible = true;
		double gas_constant = 287.1;
		double inlet_temp = 350;

		double upstream_pressure = 1000000;
		double downstream_pressure = 500000;
		double pipe_inner_dia = 0.03568248232;
		double flow_rate = 1.043534;
		double isentropic_exponent = 1.4;
		OrificeCalculation::Tapping_Option tapping_option = OrificeCalculation::Tapping_Option::Flange_Tapping;

		// Output values
		int er_code = 0;
		double orifice_dia = 0;
		OrificeCalculation test_calculation(viscosity, density,error_tolerance, compressible, 
																		gas_constant, inlet_temp);
		er_code = test_calculation.GetOrificeDiameter(upstream_pressure, downstream_pressure, flow_rate,
														  pipe_inner_dia, isentropic_exponent,
														  tapping_option, 
														  orifice_dia);

		//_ASSERT(er_code == 0 && fabs(flow_rate - 4.17) < 1e-2); // flow_rate should be 4.17
	}
    return 0;
}