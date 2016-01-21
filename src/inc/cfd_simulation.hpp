#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>

#include "defines.hpp"
#include "simulation_data.hpp"

class CFD_simulation {
	
	// class containing all the data necessary for simulating
	// a fluid by solving the 2D Navier-Stokes Equations
	
public:

	// constructor getting a pointer to the data structure for data access
	CFD_simulation(simulation_data* _dat);

	// compute one complete time-step
	void compute_step();

	//default destructor
	~CFD_simulation() = default;

private:

	simulation_data* dat;

	// computing derivatives (no index checking made for the following
	// procedures, (i,j) are supposed to be inside the right boundaries)
	
	double dux_dx(const int& i, const int& j) const;
	double duy_dy(const int& i, const int& j) const;
	
	double dux2_dx(const int& i, const int& j) const;
	double duy2_dy(const int& i, const int& j) const;
	
	double duxuy_dx(const int& i, const int& j) const;
	double duxuy_dy(const int& i, const int& j) const;
	
	double d2ux_dx2(const int& i, const int& j) const;
	double d2ux_dy2(const int& i, const int& j) const;

	double d2uy_dx2(const int& i, const int& j) const;
	double d2uy_dy2(const int& i, const int& j) const;
	
	double dp_dx(const int& i, const int& j) const;
	double dp_dy(const int& i, const int& j) const;

	// check if delta_t is still in the correct range
	void adapt_delta_t();

	// update the domain boundary conditions
	void update_domain();

	// compute F and G field
	void compute_F_G();

	// compute the righthand side of the possion equation
	void compute_Possion_RHS();

	// solve the Poisson equation using the relaxed Jacobi method
	void solve_Poisson_Jacobi();

	// compute u at the next timestep by applying the pressure correction
	void compute_u_next();
};

#endif //SIMULATION_H
