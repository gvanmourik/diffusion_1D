#ifndef CRANK_NICOLSON_H_
#define CRANK_NICOLSON_H_

#include <cmath>

// Forward time 
class Crank_Nicolson {

private:
	double K;		
	double tau;		// time step
	double h;		// spatial step
	double alpha;	// method is stable if alpha < 1/2

public:
	Crank_Nicolson(): K(1.0), tau(0.05), h(0.5) 
	{
		calculate_alpha();
	}
	Crank_Nicolson(double user_K, double user_tau, double user_h)
	{
		calculate_alpha();
	}
	~Crank_Nicolson();


	void calculate_alpha()
	{
		alpha = K * tau / pow(h,2);
	}

	
};




#endif /* CRANK_NICOLSON_H_ */