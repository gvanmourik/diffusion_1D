#include <iostream>
#include <string>
#include <FTCS.h>
#include <Crank_Nicolson.h>

void print_vector_double(std::vector<double> &src, std::string v_name);


int main(int argc, char const *argv[])
{
	std::vector<double> U_results;
	FTCS diffusion;
	diffusion.print_A_();
	diffusion.ftcs_explicit_method();
	U_results = diffusion.get_U_();

	print_vector_double(U_results, "U_results");


	return 0;
}

void print_vector_double(std::vector<double> &src, std::string v_name)
{
	printf( "\n%s Vector:\n", v_name.c_str() );
	for (int i = 0; i < src.size(); i++)
	{
		printf( "%s[%d] = %f\n", v_name.c_str(), i, src[i] );
	}
}
