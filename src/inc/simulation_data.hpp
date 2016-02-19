#ifndef SIMDATA_H
#define SIMDATA_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <mpi.h>

#include "defines.hpp"

class simulation_data {
	
	// class containing all the data necessary for doing a CFD simulation
    friend class CFD_simulation;
	
public:

	// standard constructor, filling the data with default values
	simulation_data(const std::array<double, 2> _length = {0.0, 0.0}, const std::array<unsigned int, 2> _nmax = {0, 0}, int _rank = 0, int _size = 1);
	
	// set initial conditions
	void initial_conditions(const double& u_x_init = 0.0, const double& u_y_init = 0.0, const double& p_init = 0.0);
	
	// set boundary conditions
	bool boundary_conditions(std::vector<std::vector<int>> obstacle_matrix);
	
	// write VTK file for visualisation purposes
	void export_vtk(const std::string filename) const;

	// default destructor
	~simulation_data() = default;

public:

	// --- Geometry data -------------------------------------------------
	std::array<double, 2> length;			// 2D domain size
	std::array<unsigned, 2> nmax;		    // number of interior cells
	std::array<double, 2> delta;			// grid step size

	// --- Time-stepping data --------------------------------------------
	double t;								// current time value
	double t_end;							// final time value
	double delta_t;							// time step size
	bool delta_t_is_maximal_value;			// the given delta_t is the maximal possible timestep?
	double tau;								// safety factor for time step size control

	// --- Pressure-iteration data ---------------------------------------
	unsigned itermax;					    // maximal number of pressure iterations in one time step
	double eps;								// stopping tolerance for the pressure iteration
	double omega;							// relaxation factor for the pressure iteration
	double gamma;							// upwind difference factor

	// --- Problem-dependent quantities ----------------------------------
	double Re;								// Reynolds number
	std::array<double, 2> g;				// body forces (e.g. gravity)
	std::array<int, 4> domain_boundary;		// boundary condition types for the 4 domain boundaries
	std::array<double, 2> u_max_abs;		// absolute maximum velocties appearing during one timestep
	std::array<double, 4> u_inflow;			// set inflow boundary conditions (only orthogonal)
    
    // --- Parallel computation properties ---------------------------------
    int rank;                               // id of processor
    int size;                               // number of processors
    MPI_Status status;                      // MPI Status

private:

	// --- CFD-data arrays -----------------------------------------------
	std::vector<std::vector<double>> u_x;	// velocity field in x-direction
	std::vector<std::vector<double>> u_y;	// velocity field in y-direction
	std::vector<std::vector<double>> p;    	// pressure field
    std::vector<std::vector<double>> p_new;	// updated pressure field
	std::vector<std::vector<double>> rhs;	// righthand side for pressure iteration
	std::vector<std::vector<double>> F;		// F array (abbreviation)
	std::vector<std::vector<double>> G;		// G array (abbreviation)
	std::vector<std::vector<int>> flag;		// flag field for storing boundary obstacles in the fluid field
};

#endif //SIMDATA_H
