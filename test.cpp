#include "OrificeCalculation.h"
#include <iostream>

// basic file operations
#include <iostream>
#include <fstream>
using namespace std;

int TestHardCodedCases();
int TestFileInputCases();
int TestFile();
void ExplainErrorCode(int err_code);

int main()
{
	TestHardCodedCases();
	//TestFileInputCases();

    return 0;
}

int TestFileInputCases()
{
	cout << "-----------------------------------\n";
	cout << "Testing utility for calc functions\n";
	cout << "-----------------------------------\n";
	char option = 'i';
	while (option != 'x')
	{
		cout << "\nType i to specify input file, h for help, x for exit\n";
		cin >> option;
		if (option == 'i')
		{
			TestFile();
		}
		else if (option == 'h')
		{

		}
		else if (option == 'x')
			return 0;
	}
	return 0;
}

int TestFile()
{
	cout << "\nType name of the text file, i.e. 'input.txt' and hit enter\n";

	char line[100];
	char calc_type = '0';
	char filename[100];
	cin >> filename;
	ifstream myfile (filename);
	double viscosity = -1;
	double density = -1;
	double error_tolerance = -1;
	bool compressible = true;
	double gas_constant =-1;
	double inlet_temp = -1;

	double upstream_pressure = -1;
	double downstream_pressure = -1;
	double pipe_inner_dia = -1;
	double orifice_dia = -1;
	double isentropic_exponent = -1;
	int tapping_option = -1;

	double flow_rate = -1;
	int err_code = 0;

	// Output values
	int er_code = 0;
	int no_of_inputs = 0;
	if (myfile.is_open())
	{
		while ( myfile.getline(line, 100) )
		{
			std::string str1 = line;
			int k = -1;
			k = (int) str1.find('*');
			if (k != -1)
				continue; // comment line
			else if (str1.compare("calculation_type") == 0)
			{
				myfile.getline(line, 100);
				calc_type = line[0];
			}
			else if (str1.compare("viscosity") == 0)
			{
				myfile.getline(line, 100);
				viscosity = ::atof(line);
			}
			else if (str1.compare("density") == 0)
			{
				myfile.getline(line, 100);
				density = ::atof(line);
			}
			else if (str1.compare("error_tolerance") == 0)
			{
				myfile.getline(line, 100);
				error_tolerance = ::atof(line);
			}
			else if (str1.compare("compressible") == 0)
			{
				myfile.getline(line, 100);
				compressible = (line[0] == 't')? true:false;
			}
			else if (str1.compare("gas_constant") == 0)
			{
				myfile.getline(line, 100);
				gas_constant = ::atof(line);
			}
			else if (str1.compare("inlet_temp") == 0)
			{
				myfile.getline(line, 100);
				inlet_temp = ::atof(line);
			}
			else if (str1.compare("upstream_pressure") == 0)
			{
				myfile.getline(line, 100);
				upstream_pressure = ::atof(line);
			}
			else if (str1.compare("downstream_pressure") == 0)
			{
				myfile.getline(line, 100);
				downstream_pressure = ::atof(line);
			}
			else if (str1.compare("pipe_inner_dia") == 0)
			{
				myfile.getline(line, 100);
				pipe_inner_dia = ::atof(line);
			}
			else if (str1.compare("orifice_dia") == 0)
			{
				myfile.getline(line, 100);
				orifice_dia = ::atof(line);
			}
			else if (str1.compare("flow_rate") == 0)
			{
				myfile.getline(line, 100);
				flow_rate = ::atof(line);
			}
			else if (str1.compare("isentropic_exponent") == 0)
			{
				myfile.getline(line, 100);
				isentropic_exponent = ::atof(line);
			}
			else if (str1.compare("tapping_option") == 0)
			{
				myfile.getline(line, 100);
				std::string str2 = line;
				if (!str2.compare("na"))
					tapping_option = OrificeCalculation::Not_Applicable;
				else if (!str2.compare("corner"))
					tapping_option = OrificeCalculation::Corner_Tapping;
				else if (!str2.compare("d_and_d_by_2"))
					tapping_option = OrificeCalculation::D_And_D_by_2_Tapping;
				else if (!str2.compare("flange"))
					tapping_option = OrificeCalculation::Flange_Tapping;
			}
		}
		myfile.close();
		if (viscosity < 0 || density < 0 || error_tolerance < 0 || gas_constant < 0 || inlet_temp < 0
			|| upstream_pressure < 0 || downstream_pressure < 0 || pipe_inner_dia < 0
			|| (orifice_dia < 0 && flow_rate < 0) || isentropic_exponent < 0 || calc_type == '0'
			|| tapping_option < 0)
			cout << "Input data not complete\n";
		else
			cout << "Input data read successfully\n";

		OrificeCalculation test_calculation(viscosity, density,error_tolerance, 
												compressible, gas_constant, inlet_temp);
		if (calc_type == 'f')
		{
			er_code = test_calculation.GetOrificeMassFlowRate(upstream_pressure, downstream_pressure, 
														  pipe_inner_dia, orifice_dia, isentropic_exponent,
														  tapping_option, 
														  flow_rate);
			if (err_code == 0)
				cout << "Flow Rate Calculated = " << flow_rate << "\n";
			else
				ExplainErrorCode(err_code);
		}
		else if (calc_type == 'd')
		{
			er_code = test_calculation.GetOrificeDiameter( upstream_pressure,downstream_pressure,flow_rate,
												 pipe_inner_dia, isentropic_exponent, tapping_option, orifice_dia);
			if (err_code == 0)
				cout << "Orifice Dia Calculated = " << orifice_dia << "\n";
			else
				ExplainErrorCode(err_code);
		}
		cout << "Test Case Executed Successfully" << "\n";

	}

	else cout << "Unable to open file"; 
	return 0;
}

int TestHardCodedCases()
{
	// For testing and validation, the following website was used
	// http://www.pipeflowcalculations.net/orifice.xhtml
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
		double flow_rate = 0.2498;
		double isentropic_exponent = 1.3;
		OrificeCalculation::Tapping_Option tapping_option = OrificeCalculation::Tapping_Option::Corner_Tapping;

		// Output values
		int er_code = 0;
		double diameter = 0;

		OrificeCalculation test_calculation(viscosity, density,error_tolerance, compressible, 
																		gas_constant, inlet_temp);	
		er_code = test_calculation.GetOrificeDiameter( upstream_pressure,downstream_pressure,flow_rate,
												 pipe_inner_dia, isentropic_exponent, tapping_option, diameter);
		_ASSERT(er_code == 0 &&  fabs(diameter - 0.02489) < 1e-4); // diameter shold be 0.02480	
	}

	// Test case 3 : Orifice flowrate, incompressible
	if (false)
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
	if (false)
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
														  orifice_dia);
		_ASSERT(er_code == 0 && fabs(orifice_dia - 0.025) < 1e-2); // flow_rate should be 0.025
	}

	// Test case 5 : Orifice flowrate, tapping option - flange
	if(true)
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

		_ASSERT(er_code == 0 && fabs(flow_rate - 0.2443) < 1e-2); 
	}

	// Test case 6 : Choked flow, compressible, calculate mass flow rate
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

		_ASSERT(er_code == 0 && fabs(flow_rate - 1.0434) < 1e-2); 
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

		_ASSERT(er_code == 0 && fabs(orifice_dia - 0.0248) < 1e-2); 
	}
	return 0;
}

void ExplainErrorCode(int err_code)
{
	switch (err_code)
	{
		case ERCODE_INVALID_INPUT:
			{
				cout << "Inputs are invalid \n";
			}
			break;
		case ERCODE_MAX_ITERATIONS_EXCEED:
			{
				cout << "Maximum iterations exceeded in calculation \n";
			}
			break;
		case ERCODE_ERROR_DIVERGING:
			{
				cout << "Error in iterative calculation is diverging\n";
			}
			break;		
		case ERCODE_PRESSURE_DROP_NEGATIVE:
			{
				cout << "Pressuredrop is negative \n";
			}
			break;		
		case ERCODE_MACH_NUMBER_NEGATIVE:
			{
				cout << "Mach number calculated is negative \n";
			}
			break;		
		case ERCODE_REYNOLDS_NUMBER_NEGATIVE:
			{
				cout << "reynolds number calculated is negative \n";
			}
			break;		
		case ERCODE_DISCHARGE_COEFFICIENT_NEGATIVE:
			{
				cout << "Discharge coefficient calculated is negative \n";
			}
			break;
	}
}