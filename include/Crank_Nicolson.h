#ifndef CRANK_NICOLSON_H_
#define CRANK_NICOLSON_H_

#include <cmath>
#include <string>
#include <vector>


class Crank_Nicolson {

private:
	int N;									// defines dimensions of A and x
	double K;		
	double tau;								// time step
	double h;								// spatial step
	double alpha;							// method is stable if alpha < 1/2
	double d;								// Main diagonal elements value
	double time_steps;
	std::vector<double> U;					// Vector of temperatures across the bar
	std::vector<double> U_star;
	std::vector<double> c_prime, U_prime;	// Used in the Thomas Algorithm
	std::vector<double> a, b, c;			// Sub, main, and super-diagonal components
	                            			// 	|	\	\	\		 |
	                            			// 	|	 \	 \	 \		 |
	                            			// 	|	  a	  b	  c		 |
	                            			// 	|	   \   \   \	 |
	                            			// 	|	 	\	\	\	 |

public:
	Crank_Nicolson(): alpha(0.499), h(0.1), K(1.0) 
	{
		set_tau();
		set_d();
		set_element_dimensions();
		set_diagonal_values(); // for a, b, and c
		set_initial_U_();
	}
	Crank_Nicolson(double user_alpha, double user_h, double user_K)
	{
		alpha = user_alpha;
		h = user_h;
		K = user_K;

		set_tau();
		set_d();
		set_element_dimensions();
		set_diagonal_values(); // for a, b, and c
		set_initial_U_();
	}
	~Crank_Nicolson(){}


	// Explicit FTCS Method
	void crank_nicolson_method()
	{

		printf("-->Time Steps = %f\n", time_steps);
		for (int i = 0; i < time_steps; i++)
		{
			set_U_star_();
			U = U_next_();
			
			print_vector(U, "U_next");
			printf("\n");
		}
	}

	// Thomas algorithm to find the next U
	std::vector<double> U_next_()
	{
		std::vector<double> U_next(N);

		double m;
		for (int i = 1; i < N-1; i++)
		{
			m = a[i-1] / b[i-1];
			b[i] -= (m * c[i-1]);
			U_star[i] -= (m * U_star[i-1]);
		}

		// Solve for the last U_next value
		U_next[N-1] = U_star[N-1] / b[N-1];

		// printf("U_star[N-1] = %f\n", U_star[N-1]);
		// printf("b[N-1] = %f\n", b[N-1]);

		// Solve for the rest using back substitution
		for (int i = N-2; i >= 0; i--)
		{
			U_next[i] = ( U_star[i] - (c[i]*U_next[i+1]) ) / b[i];
		}

		return U_next;
	}

	void set_U_star_()
	{
		int n = N - 1;

		U_star[0] = U[1] * alpha/2;
		U_star[n] = U[n-1] * alpha/2;

		for (int i = 1; i < n; i++)
		{
			U_star[i] = (U[i-1] * alpha/2.0) + (U[i] * (1.0-alpha)) + (U[i+1] * alpha/2.0);

			// printf("U[%d] * alpha/2.0 = %f\n", i-1, U[i-1] * alpha/2.0);
			// printf("U[%d] * (1-alpha) = %f\n", i, U[i] * (1.0-alpha));
			// printf("U[%d] * alpha/2.0 = %f\n\n", i+1, U[i+1] * alpha/2.0);
		}

		print_vector(U_star, "U_star");
	}

	// Initial U Vector
	void set_initial_U_()
	{
		// Set spatial boundary conditions
		U[0] = 0.0;
		U[N-1] = 0.0;

		printf("\nN = %d\n\n", N);
		for (int i = 1; i < N-1; i++)
		{
			U[i] = U_x_0_(i * h);
			// printf("i * h = %f\n", i * h);
		}

		print_vector(U, "Initial U_vector");
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
		time_steps = (1/tau) + 1;
	}

	void set_d()
	{
		d = 1.0 + alpha;
	}

	void set_diagonal_values()
	{
		a.resize(N);
		b.resize(N);
		c.resize(N);

		for (int i = 0; i < N; i++)
		{
			a[i] = -alpha/2;
			b[i] = d;
			c[i] = -alpha/2;
		}
		
		print_vector(a, "A");
		print_vector(b, "B");
		print_vector(c, "C");
	}

	// void set_prime_vectors()
	// {
	// 	c_prime.resize(N);
	// 	U_prime.resize(N);

	// 	c_prime[0] = c[0] / b[0];
	// 	U_prime[0] = d[0] / b[0];
	// 	double beta;

	// 	for (int i = 1; i < N-1; i++)
	// 	{
	// 		beta = 1 / ( b[i] - c_prime[i-1] * a[i] );

	// 		c_prime[i] = beta * c[i];
	// 		U_prime[i] = beta * ( U_prime[i] - U_prime[i-1] * a[i] );
	// 	}
	// }

	void set_element_dimensions()
	{
		// Set N
		N = (1/h) + 1;
		// Resize U according to N
		U.resize(N);
		U_star.resize(N);
	}

	void set_alpha()
	{
		alpha = K * tau / pow(h,2);
	}

	std::vector<double> get_U_()
	{
		return U;
	}	

	void print_vector(std::vector<double> &src, std::string v_name)
	{
		printf( "%s:\n", v_name.c_str());
		printf("[   ");
		for (int i = 0; i < src.size(); i++)
		{
			printf("%f   ", src[i]);
		}
		printf("]\n");
	}


};


#endif /* CRANK_NICOLSON_H_ */