#ifndef CRANK_NICOLSON_H_
#define CRANK_NICOLSON_H_

#include <cmath>
#include <string>
#include <vector>
#include <fstream>


class Crank_Nicolson {

private:
	int N;									// defines dimensions of A and x
	int time_steps;
	double K;		
	double tau;								// time step
	double h;								// spatial step
	double alpha;							// method is stable if alpha < 1/2
	double d;								// Main diagonal elements value
	std::vector<double> U;					// Vector of temperatures across the bar
	std::vector<double> U_star;				// U^* vector defined in the report
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
		//set_prime_vectors();
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
		//set_prime_vectors();
		set_initial_U_();
	}
	~Crank_Nicolson(){}


	// Explicit FTCS Method
	void crank_nicolson_method()
	{

		// printf("\n-->Time Steps = %d\n\n", time_steps);
		for (int i = 0; i < time_steps; i++)
		{
			set_U_star_();
			U = U_next_();
		}
	}

	// Thomas algorithm to find the next U
	std::vector<double> U_next_()
	{
		int n = N - 1;
		set_prime_vectors();
		std::vector<double> U_next(N);

		double U_prime_n_ = ( U_star[n] - U_prime[n-1] * a[n] ) /
							         ( b[n] - c_prime[n-1] * a[n] );
		
		// Initial U_next value, followed by backwards substitution.
		U_next[n] = U_prime_n_;
		for (int i = n-1; i >= 0; i--)
		{
			U_next[i] = U_prime[i] - c_prime[i] * U_next[i+1];
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
		}
	}

	// Initial U Vector
	void set_initial_U_()
	{
		// Set spatial boundary conditions
		U[0] = 0.0;
		U[N-1] = 0.0;

		// printf("\nN = %d\n\n", N);
		for (int i = 1; i < N-1; i++)
		{
			U[i] = U_x_0_(i * h);
		}
	}

	double U_x_0_(double x)
	{
		if (0.0 < x && x <= 0.5)
		{
			return 2.0 * x;
		}
		if (0.5 < x && x <= 1.0)
		{
			return 2.0 * (1.0-x);
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
			a[i] = -alpha/2.0;
			b[i] = d;
			c[i] = -alpha/2.0;
		}
	}

	void set_prime_vectors()
	{
		c_prime.resize(N);
		U_prime.resize(N);

		c_prime[0] = c[0] / b[0];
		U_prime[0] = U_star[0] / b[0];
		double beta;

		for (int i = 1; i < N; i++)
		{
			beta = 1.0 / ( b[i] - c_prime[i-1] * a[i] );

			c_prime[i] = beta * c[i];
			U_prime[i] = beta * ( U_star[i] - U_prime[i-1] * a[i] );
		}
	}

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

	void normalize_vector(std::vector<double> &src)
	{
		double max = *std::max_element(src.begin(), src.end());

		for (int i = 0; i < N; i++)
		{
			src[i] = src[i] / max;
		}
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

	bool export_U_vs_t_(std::vector<double> &x, std::string fname)
	{
		if ( N != x.size() )
		{
			printf("Please resize %s to N elements and try again...\n", fname.c_str());
			return false;
		}

		std::ofstream file("data/" + fname);	
		normalize_vector(U);
		if ( file.is_open() )
	    {
	       	for (int i = 0; i < N; i++)
	       	{
	        	file << x[i] << " " << U[i] << "\n";
	        }
	        return true;
	    }
	    else
	    {
	    	perror("ERROR Unable to export to CSV");
	    }

	    return false;
	}


};

#endif /* CRANK_NICOLSON_H_ */