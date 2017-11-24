#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>


const int SUDO_INF = 10000; 	//max int range


double U_x_t_(double t);
std::vector<double> U_x_t_values(int steps, std::vector<double> &t);
std::vector<double> generate_t_(int steps, double t_min, double t_max);
void export_xy_to_file(std::vector<double> x, std::vector<double> y, std::string fname);


int main(int argc, char* argv[])
{
	int steps = 100;
	std::vector<double> U, t;
	t = generate_t_(steps, 0.0, 1.0);

	U = U_x_t_values(steps, t);
	printf("exporting U vs t data... {data/U_vs_t_analytical.data}\n");
	export_xy_to_file(t, U, "U_vs_t_analytical.data");


	return 0;
}

double U_x_t_(double t)
{
	double sum = 0.0;
	double pi = M_PI; 
	double exp_arg;

	for (int n = 1; n < SUDO_INF; n++)
	{
		exp_arg = pow(n*pi,2) * t;
		// printf("exp_arg[%d] = %f\n", n, exp_arg);

		sum += (1/pow(n,2)) * sin(n*pi/2) * exp(-exp_arg);
		// printf("sum = %f\n", sum);
	}

	return sum * 8.0 / pow(pi,2);
}

std::vector<double> U_x_t_values(int steps, std::vector<double> &t)
{
	std::vector<double> U(steps);

	for (int i = 0; i < steps; i++)
	{
		U[i] = U_x_t_(t[i]);
		// printf("U[%d] = %f\n", i, U[i]);
	}

	return U;
}

std::vector<double> generate_t_(int steps, double t_min, double t_max)
{
	double dt = (t_max - t_min) / (steps-1);
	std::vector<double> t(steps);

	for (int i = 0; i < steps; i++)
	{
		t[i] = t_min + (i * dt);
	}

	return t;
}

void export_xy_to_file(std::vector<double> x, std::vector<double> y, std::string fname)
{
	std::ofstream file("data/" + fname);	

	if ( file.is_open() )
    {
  	 	int length = y.size();
       	for (int i = 0; i < length; i++)
       	{
        	file << x[i] << " " << y[i] << "\n";
        }
    }
    else
    {
    	perror("ERROR Unable to export to CSV");
    }
}
