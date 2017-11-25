#include <iostream>
#include <string>
#include <FTCS.h>
#include <Crank_Nicolson.h>

void print_vector_double(std::vector<double> &src, std::string v_name);


int main(int argc, char const *argv[])
{
	std::vector<double> U_results;


	// FTCS diffusion_stable(0.499, 0.1, 1.0);
	// diffusion_stable.ftcs_explicit_method();
	// U_results = diffusion_stable.get_U_();
	// print_vector_double(U_results, "Explicit Stable Results");

	// FTCS diffusion_unstable(0.502, 0.1, 1.0);
	// diffusion_unstable.ftcs_explicit_method();
	// U_results = diffusion_unstable.get_U_();
	// print_vector_double(U_results, "Explicit Unstable Results");

	Crank_Nicolson diffusion_stable_cn(0.499, 0.1, 1.0);
	diffusion_stable_cn.crank_nicolson_method();
	U_results = diffusion_stable_cn.get_U_();
	print_vector_double(U_results, "CN Stable Results");

	// Crank_Nicolson diffusion_unstable_cn(0.502, 0.1, 1.0);
	// diffusion_unstable_cn.crank_nicolson_method();
	// U_results = diffusion_unstable_cn.get_U_();
	// print_vector_double(U_results, "CN Unstable Results");


	return 0;
}

void print_vector_double(std::vector<double> &src, std::string v_name)
{
	printf( "\n%s:\n", v_name.c_str() );
	for (int i = 0; i < src.size(); i++)
	{
		printf( "%s[%d] = %f\n", v_name.c_str(), i, src[i] );
	}
}
