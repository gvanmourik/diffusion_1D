#include <iostream>
#include <string>
#include <FTCS.h>
#include <Crank_Nicolson.h>

void print_vector_double(std::vector<double> &src, std::string v_name);
std::vector<double> generate_t_(int size, double t_min, double t_max);

int main(int argc, char const *argv[])
{
	int spatial_steps;
	double step_size = 0.01;
	std::vector<double> U_numerical, U_analytical, t;

	FTCS diffusion_stable(0.499, step_size, 1.0);
	U_numerical = diffusion_stable.get_U_();
	spatial_steps = U_numerical.size();
	t = generate_t_(spatial_steps, 0, 1);
	diffusion_stable.export_U_vs_t_(t, U_numerical, "initial_distribution.data", false);
	
	diffusion_stable.ftcs_explicit_method();
	U_numerical = diffusion_stable.get_U_();
	t = generate_t_(spatial_steps, 0, 1);
	U_analytical = diffusion_stable.U_analytical(spatial_steps, t);
	diffusion_stable.export_U_vs_t_(t, U_numerical, "FTCS_stable.data", true);
	diffusion_stable.export_U_vs_t_(t, U_analytical, "Analytical.data", true);
	
	FTCS diffusion_unstable(0.502, step_size, 1.0);
	diffusion_unstable.ftcs_explicit_method();
	U_numerical = diffusion_unstable.get_U_();
	// diffusion_unstable.fill_w_avg_vector_value(U_numerical);
	t = generate_t_(U_numerical.size(), 0, 1);
	diffusion_unstable.export_U_vs_t_(t, U_numerical, "FTCS_unstable.data", true);

	Crank_Nicolson diffusion_CN(0.499, step_size, 1.0);
	diffusion_CN.crank_nicolson_method();
	U_numerical = diffusion_CN.get_U_();
	t = generate_t_(U_numerical.size(), 0, 1);
	diffusion_CN.export_U_vs_t_(t, "CN.data");


	return 0;
}

std::vector<double> generate_t_(int size, double t_min, double t_max)
{
	double tau = (t_max - t_min) / (size - 1);
	std::vector<double> t(size);

	for (int i = 0; i < size; i++)
	{
		t[i] = t_min + (i * tau);
	}

	return t;
}

void print_vector_double(std::vector<double> &src, std::string v_name)
{
	printf( "\n%s:\n", v_name.c_str() );
	for (int i = 0; i < src.size(); i++)
	{
		printf( "%s[%d] = %f\n", v_name.c_str(), i, src[i] );
	}
}
