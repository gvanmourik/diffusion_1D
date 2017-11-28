#ifndef FTCS_H_
#define FTCS_H_

#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>


// Forward-time centered-space 
class FTCS {

private:
	int N;									// defines dimensions of A and x
	double K;		
	double tau;								// time step
	double h;								// spatial step
	double alpha;							// method is stable if alpha < 1/2
	double d;								// Main diagonal elements value
	std::vector<double> U;					// Vector of temperatures across the bar

public:
	FTCS(): alpha(0.499), h(0.1), K(1.0) 
	{
		set_tau();
		set_d();
		set_element_dimensions();
		construct_initial_U_();
	}
	FTCS(double user_alpha, double user_h, double user_K)
	{
		alpha = user_alpha;
		h = user_h;
		K = user_K;

		set_tau();
		set_d();
		set_element_dimensions();
		construct_initial_U_();
	}
	~FTCS(){}


	// Explicit FTCS Method
	void ftcs_explicit_method()
	{
		int time_steps = (1/tau) + 1;

		for (int i = 0; i < time_steps; i++)
		{
			U = U_next_();
		}
	}

	std::vector<double> U_next_()
	{
		int n = N - 1;
		std::vector<double> U_next(N);

		// Set values at BC's
		U_next[0] = U[0] * 1.0;
		U_next[n] = U[n] * 1.0;

		for (int i = 1; i < n; i++)
		{
			U_next[i] = (U[i-1] * alpha) + (U[i] * d) + (U[i+1] * alpha);
		}

		return U_next;
	}

	// Initial U Vector
	void construct_initial_U_()
	{
		// Set spatial boundary conditions
		U[0] = 0.0;
		U[N-1] = 0.0;

		for (int i = 1; i < N-1; i++)
		{
			U[i] = U_x_0_(i * h);
		}
	}

	double U_x_0_(double x)
	{
		if (0.0 < x && x <= 0.5)
		{
			return 2 * x;
		}
		if (0.5 < x && x <= 1.0)
		{
			return 2 * (1-x);
		}

		return 0.0;
	}


	// Other
	void set_tau()
	{
		tau = alpha * pow(h,2) / K;
	}

	void set_d()
	{
		d = 1.0 - (2.0 * alpha);
	}

	void set_element_dimensions()
	{
		// Set N
		N = (1/h) + 1;

		// Resize U according to N
		U.resize(N);
	}

	std::vector<double> get_U_()
	{
		return U;
	}

	void normalize_vector(std::vector<double> &src)
	{
		double max = *std::max_element(src.begin(), src.end());

		for (int i = 0; i < N; i++)
		{
			src[i] = src[i] / max;
		}
	}

	void fill_w_avg_vector_value(std::vector<double> &src)
	{
		normalize_vector(src);

		// Find the average value of the vector
		int elements = src.size();
		double avg_value, sum = 0.0;
		for (int i = 0; i < elements; i++)
		{
			sum += src[i];
		}
		avg_value = sum / (double)elements;

		for (int i = 0; i < elements; i++)
		{
			src[i] = avg_value;
		}
	}	

	void print_U_()
	{
		printf( "U Vector:\n");
		printf("[   ");
		for (int i = 0; i < U.size(); i++)
		{
			printf("%f   ", U[i]);
		}
		printf("]\n");
	}

	bool export_U_vs_t_(std::vector<double> &x, std::vector<double> y, std::string fname, bool norm)
	{
		if ( N != x.size() )
		{
			printf("Please resize %s to N elements and try again...\n", fname.c_str());
			return false;
		}

		if (norm)
		{
			normalize_vector(y);
		}	

		std::ofstream file("data/" + fname);
		if ( file.is_open() )
	    {
	       	for (int i = 0; i < N; i++)
	       	{
	        	file << x[i] << " " << y[i] << "\n";
	        }
	        return true;
	    }
	    else
	    {
	    	perror("ERROR Unable to export to CSV");
	    }

	    return false;
	}

	// Analytical Solution
	std::vector<double> U_analytical(int steps, std::vector<double> &t)
	{
		std::vector<double> U_analytical(steps);

		for (int i = 0; i < N; i++)
		{
			U_analytical[i] = U_x_t_(i*h, t);
			// printf("i*h = %f\n", i*h);
		}

		return U_analytical;
	}

	double U_x_t_(double x, std::vector<double> &t)
	{
		double sum = 0.0;
		double pi = M_PI; 
		double exp_arg;

		for (int n = 1; n < N; n++)
		{
			exp_arg = -pow(n*pi,2) * t[n];
			// printf("exp_arg[%d] = %f\n", n, exp_arg);

			sum += (1.0/pow(n,2)) * sin(n*pi/2) * sin(n*pi*x) * exp(exp_arg);
			// printf("sum = %f\n", sum);
		}

		return sum * 8.0 / pow(pi,2);
	}


};

#endif /* FTCS_H_ */