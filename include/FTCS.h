#ifndef FTCS_H_
#define FTCS_H_

#include <cmath>
#include <vector>


extern "C" {
     void dgemv_(char *trans, int *M, int *N, double *my_alpha, 
     	double *A, int *lda, double *x, int *incx, double *beta, 
     	double *y, int *incy, int *error);
}


// Forward time 
class FTCS {

private:
	int N;									// defines dimensions of A and x
	double K;		
	double tau;								// time step
	double h;								// spatial step
	double alpha;							// method is stable if alpha < 1/2
	std::vector<double> U;					// Vector of temperatures across the bar
	std::vector< std::vector<double> > A;	// A matrix defined in the report

public:
	FTCS(): alpha(0.499), h(0.01), K(1.0) 
	{
		set_tau();
		set_element_dimensions();
		construct_A_();
		construct_initial_U_();
	}
	FTCS(double user_alpha, double user_h, double user_K)
	{
		alpha = user_alpha;
		h = user_h;
		K = user_K;

		set_tau();
		set_element_dimensions();
		construct_A_();
		construct_initial_U_();
	}
	~FTCS(){}


	// Explicit FTCS Method
	void ftcs_explicit_method()
	{
		int iterations = 1.0 / tau;

		printf("\n");
		for (int i = 0; i < 15; i++)
		{
			U = U_next();
			print_U_();
		}
	}

	// ISSUES with BLAS and LAPACK. Only works with element size of 2^n.
	//
	// std::vector<double> U_next()
	// {
	// 	// Values for dtrmv_
	// 	char trans;
	// 	double my_alpha, beta;
	// 	int M, lda, incx, incy, error;
	// 	std::vector<double> U_next(N);

	// 	trans = 't';
	// 	M = N;
	// 	my_alpha = 1.0;
	// 	beta = 0.0;
	// 	lda = M;
	// 	incx = 1;
	// 	incy = 1;

	// 	// printf("N = %d, M = %d\n", N, M );
		

	// 	dgemv_(&trans, &M, &N, &my_alpha, &A[0][0], &lda, &U[0], &incx, &beta, &U_next[0], &incy, &error);
	// 	// dgemv_(&trans, &M, &N, &alpha, &A[0][0], &lda, &x[0], &incx, &beta, &y[0], &incy, &error);

	// 	return U_next;
	// }


	// A Matrix
	void construct_A_()
	{	
		boundary_conditions_A_();

		double diag_value = 1 - 2*alpha;
		double sub_diag_value = alpha;

		for (int i = 1; i < N-1; i++)
		{
			A[i][i] = diag_value;
			A[i][i-1] = sub_diag_value;
			A[i][i+1] = sub_diag_value;
		}

	}

	void boundary_conditions_A_()
	{
		int n = N-1;

		A[0][0] = 1.0;
		A[n][n] = 1.0;

		for (int i = 1; i < N; i++)
		{
			A[0][i] = 0.0;
			A[n][n-i] = 0.0;
		}
	}


	// Initial T Vector
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

	void set_element_dimensions()
	{
		// Set N
		N = (1/tau) + 1;

		// Initialize A and U according to N
		A.resize( N, std::vector<double>(N,0.0) );
		U.resize( N );
	}

	void set_alpha()
	{
		alpha = K * tau / pow(h,2);
	}

	std::vector<double> get_U_()
	{
		return U;
	}

	void print_A_()
	{
		printf("\nMatrix A:\n");
		for (int i = 0; i < N; i++)
		{
			printf("[%d] [   ", i);
			for (int j = 0; j < N; j++)
			{
				printf("%f   ", A[i][j]);
			}
			printf("]\n");
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


};




#endif /* FTCS_H_ */