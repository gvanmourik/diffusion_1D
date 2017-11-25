#ifndef FTCS_H_
#define FTCS_H_

#include <cmath>
#include <vector>


// Forward time centered space 
class FTCS {

private:
	int N;									// defines dimensions of A and x
	double K;		
	double tau;								// time step
	double h;								// spatial step
	double alpha;							// method is stable if alpha < 1/2
	double d;								// Main diagonal elements value
	std::vector<double> U;					// Vector of temperatures across the bar
	// std::vector< std::vector<double> > A;	// A matrix defined in the report

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

		printf("\n");
		for (int i = 0; i < time_steps; i++)
		{
			U = U_next_();
			//print_U_();
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

		printf("N = %d\n", N);
		for (int i = 1; i < N-1; i++)
		{
			U[i] = U_x_0_(i * h);
			// printf("i * h = %f\n", i * h);
		}

		printf("Initial U:\n");
		print_U_();
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
		// A.resize( N, std::vector<double>(N,0.0) );
		U.resize(N);
	}

	void set_alpha()
	{
		alpha = K * tau / pow(h,2);
	}

	std::vector<double> get_U_()
	{
		return U;
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


};


#endif /* FTCS_H_ */