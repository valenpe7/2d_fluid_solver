#include "inc/cfd_simulation.hpp"
#include <mpi.h>

CFD_simulation::CFD_simulation(simulation_data* _dat)
	: dat(_dat) {
}

double CFD_simulation::dux_dx(const int& i, const int& j) const {
	return (dat->u_x[i][j] - dat->u_x[i - 1][j]) / dat->delta[X];
}

double CFD_simulation::duy_dy(const int& i, const int& j) const {
	return (dat->u_y[i][j] - dat->u_y[i][j - 1]) / dat->delta[Y];
}

double CFD_simulation::dux2_dx(const int& i, const int& j) const {
	return
		1.0 / dat->delta[X] *
		(((dat->u_x[i][j] + dat->u_x[i + 1][j]) / 2.0) * ((dat->u_x[i][j] + dat->u_x[i + 1][j]) / 2.0)
			- ((dat->u_x[i - 1][j] + dat->u_x[i][j]) / 2.0) * ((dat->u_x[i - 1][j] + dat->u_x[i][j]) / 2.0))
		+ dat->gamma / dat->delta[X] *
		((fabs(dat->u_x[i][j] + dat->u_x[i + 1][j]) / 2.0) * ((dat->u_x[i][j] - dat->u_x[i + 1][j]) / 2.0)
			- (fabs(dat->u_x[i - 1][j] + dat->u_x[i][j]) / 2.0) * ((dat->u_x[i - 1][j] - dat->u_x[i][j]) / 2.0));
}

double CFD_simulation::duy2_dy(const int& i, const int& j) const {
	return
		1.0 / dat->delta[Y] *
		(((dat->u_y[i][j] + dat->u_y[i][j + 1]) / 2.0) * ((dat->u_y[i][j] + dat->u_y[i][j + 1]) / 2.0)
			- ((dat->u_y[i][j - 1] + dat->u_y[i][j]) / 2.0) * ((dat->u_y[i][j - 1] + dat->u_y[i][j]) / 2.0))
		+ dat->gamma / dat->delta[Y] *
		((fabs(dat->u_y[i][j] + dat->u_y[i][j + 1]) / 2.0) * ((dat->u_y[i][j] - dat->u_y[i][j + 1]) / 2.0)
			- (fabs(dat->u_y[i][j - 1] + dat->u_y[i][j]) / 2.0) * ((dat->u_y[i][j - 1] - dat->u_y[i][j]) / 2.0));
}

double CFD_simulation::duxuy_dx(const int& i, const int& j) const {
	return
		1.0 / dat->delta[X] *
		(((dat->u_x[i][j] + dat->u_x[i][j + 1]) / 2.0) * ((dat->u_y[i][j] + dat->u_y[i + 1][j]) / 2.0)
			- ((dat->u_x[i - 1][j] + dat->u_x[i - 1][j + 1]) / 2.0) * ((dat->u_y[i - 1][j] + dat->u_y[i][j]) / 2.0))
		+ dat->gamma / dat->delta[X] *
		((fabs(dat->u_x[i][j] + dat->u_x[i][j + 1]) / 2.0) * ((dat->u_y[i][j] - dat->u_y[i + 1][j]) / 2.0)
			- (fabs(dat->u_x[i - 1][j] + dat->u_x[i - 1][j + 1]) / 2.0) * ((dat->u_y[i - 1][j] - dat->u_y[i][j]) / 2.0));
}

double CFD_simulation::duxuy_dy(const int& i, const int& j) const {
	return
		1.0 / dat->delta[Y] *
		(((dat->u_y[i][j] + dat->u_y[i + 1][j]) / 2.0) * ((dat->u_x[i][j] + dat->u_x[i][j + 1]) / 2.0)
			- ((dat->u_y[i][j - 1] + dat->u_y[i + 1][j - 1]) / 2.0) * ((dat->u_x[i][j - 1] + dat->u_x[i][j]) / 2.0))
		+ dat->gamma / dat->delta[Y] *
		((fabs(dat->u_y[i][j] + dat->u_y[i + 1][j]) / 2.0) * ((dat->u_x[i][j] - dat->u_x[i][j + 1]) / 2.0)
			- (fabs(dat->u_y[i][j - 1] + dat->u_y[i + 1][j - 1]) / 2.0) * ((dat->u_x[i][j - 1] - dat->u_x[i][j]) / 2.0));
}

double CFD_simulation::d2ux_dx2(const int& i, const int& j) const {
	return (dat->u_x[i + 1][j] - 2.0*dat->u_x[i][j] + dat->u_x[i - 1][j]) / (dat->delta[X] * dat->delta[X]);
}

double CFD_simulation::d2ux_dy2(const int& i, const int& j) const {
	return (dat->u_x[i][j + 1] - 2.0*dat->u_x[i][j] + dat->u_x[i][j - 1]) / (dat->delta[Y] * dat->delta[Y]);
}

double CFD_simulation::d2uy_dx2(const int& i, const int& j) const {
	return (dat->u_y[i + 1][j] - 2.0*dat->u_y[i][j] + dat->u_y[i - 1][j]) / (dat->delta[X] * dat->delta[X]);
}

double CFD_simulation::d2uy_dy2(const int& i, const int& j) const {
	return (dat->u_y[i][j + 1] - 2.0*dat->u_y[i][j] + dat->u_y[i][j - 1]) / (dat->delta[Y] * dat->delta[Y]);
}

double CFD_simulation::dp_dx(const int& i, const int& j) const {
	return (dat->p[i + 1][j] - dat->p[i][j]) / dat->delta[X];
}

double CFD_simulation::dp_dy(const int& i, const int& j) const {
	return (dat->p[i][j + 1] - dat->p[i][j]) / dat->delta[Y];
}

void CFD_simulation::adapt_delta_t() {
	// compute new gamma
	dat->gamma = std::max(dat->u_max_abs[X] * dat->delta_t / dat->delta[X], dat->u_max_abs[Y] * dat->delta_t / dat->delta[Y]);
	if(dat->gamma > 1.0) {
		dat->gamma = 1.0;
	}

	// compute new delta_t
	double cond1 = dat->Re / (2.0 * (1.0 / (dat->delta[X] * dat->delta[X]) + 1.0 / (dat->delta[Y] * dat->delta[Y])));
	double cond2 = std::min(dat->delta[X] / dat->u_max_abs[X], dat->delta[Y] / dat->u_max_abs[Y]);
	if(dat->delta_t_is_maximal_value) {
		dat->delta_t = std::min(dat->delta_t, dat->tau * std::min(cond1, cond2));
	} else {
		dat->delta_t = dat->tau * std::min(cond1, cond2);
	}

	// debug info
    if(dat->rank == 0) {
        std::cout << "set delta_t to " << dat->delta_t << " and gamma to " << dat->gamma << std::endl;
    }
}

void CFD_simulation::update_domain() {
	// aliases for easy access to boundary cells
	const unsigned& imax = dat->nmax[X];
	const unsigned& jmax = dat->nmax[Y];
    
    unsigned jstart = 1;
    unsigned jend = jmax;
    
    unsigned istart = dat->rank*imax/dat->size + 1;
    unsigned iend = (dat->rank + 1)*imax/dat->size;
    
    // treat west boundary -----------------------------------------------
    if(dat->rank == 0) {
        switch(dat->domain_boundary[D_WEST]) {
            case(BC_FREE_SLIP) : {
                for(std::size_t j = jstart; j <= jend; ++j) {
                    dat->u_x[0][j] = 0;
                    dat->u_y[0][j] = dat->u_y[1][j];
                }
            } break;

            case(BC_OUTFLOW) : {
                for(std::size_t j = jstart; j <= jend; ++j) {
                    dat->u_x[0][j] = dat->u_x[1][j];
                    dat->u_y[0][j] = dat->u_y[1][j];
                }
            } break;

            case(BC_INFLOW) : {
                for(std::size_t j = jstart; j <= jend; ++j) {
                    dat->u_x[0][j] = dat->u_inflow[D_WEST];		// fixed u_in
                    dat->u_y[0][j] = dat->u_y[1][j]; 			// slip-bc for v
                }
            } break;

            case(BC_NO_SLIP) :
            default: // defaults to BC_NO_SLIP
            {
                for(std::size_t j = jstart; j <= jend; ++j) {
                    dat->u_x[0][j] = 0;
                    dat->u_y[0][j] = -1.0*dat->u_y[1][j];
                }
            } break;
        }
    }
    
    // treat east boundary -----------------------------------------------
    if(dat->rank == dat->size - 1) {
        switch(dat->domain_boundary[D_EAST]) {
            case(BC_FREE_SLIP) : {
                for(std::size_t j = jstart; j <= jend; ++j) {
                    dat->u_x[imax][j] = 0;
                    dat->u_y[imax + 1][j] = dat->u_y[imax][j];
                }
            } break;

            case(BC_OUTFLOW) : {
                for(std::size_t j = jstart; j <= jend; ++j) {
                    dat->u_x[imax][j] = dat->u_x[imax - 1][j];
                    dat->u_y[imax + 1][j] = dat->u_y[imax][j];
                }
            } break;

            case(BC_INFLOW) : {
                for(std::size_t j = jstart; j <= jend; ++j) {
                    dat->u_x[imax][j] = dat->u_inflow[D_EAST];	// fixed u_in
                    dat->u_y[imax + 1][j] = dat->u_y[imax][j];	// slip-bc for v
                }
            } break;

            case(BC_NO_SLIP) :
            default: // defaults to BC_NO_SLIP
            {
                for(std::size_t j = jstart; j <= jend; ++j) {
                    dat->u_x[imax][j] = 0;
                    dat->u_y[imax + 1][j] = -1.0*dat->u_y[imax][j];
                }
            } break;
        }
    }

	// treat north boundary ----------------------------------------------
	switch(dat->domain_boundary[D_NORTH]) {
		case(BC_FREE_SLIP) : {
			for(std::size_t i = istart; i <= iend; ++i) {
				dat->u_x[i][jmax + 1] = dat->u_x[i][jmax];
				dat->u_y[i][jmax] = 0;
			}
		} break;

		case(BC_OUTFLOW) : {
			for(std::size_t i = istart; i <= iend; ++i) {
				dat->u_x[i][jmax + 1] = dat->u_x[i][jmax];
				dat->u_y[i][jmax] = dat->u_y[i][jmax - 1];
			}
		} break;

		case(BC_INFLOW) : {
			for(std::size_t i = istart; i <= iend; ++i) {
				dat->u_x[i][jmax + 1] = dat->u_x[i][jmax];		// slip-bc for u
				dat->u_y[i][jmax] = dat->u_inflow[D_NORTH]; 	// fixed v_in
			}
		} break;

		case(BC_NO_SLIP) :
		default: // defaults to BC_NO_SLIP
		{
			for(std::size_t i = istart; i <= iend; ++i) {
				dat->u_x[i][jmax + 1] = -1.0*dat->u_x[i][jmax];
				dat->u_y[i][jmax] = 0;
			}
		} break;
	}

	// treat south boundary ----------------------------------------------
	switch(dat->domain_boundary[D_SOUTH]) {
		case(BC_FREE_SLIP) : {
			for(std::size_t i = istart; i <= iend; ++i) {
				dat->u_x[i][0] = dat->u_x[i][1];
				dat->u_y[i][0] = 0;
			}
		} break;

		case(BC_OUTFLOW) : {
			for(std::size_t i = istart; i <= iend; ++i) {
				dat->u_x[i][0] = dat->u_x[i][1];
				dat->u_y[i][0] = dat->u_y[i][1];
			}
		} break;

		case(BC_INFLOW) : {
			for(std::size_t i = istart; i <= iend; ++i) {
				dat->u_x[i][0] = dat->u_x[i][1];			// slip-bc for u
				dat->u_y[i][0] = dat->u_inflow[D_SOUTH];	// fixed v_in
			}
		} break;

		case(BC_NO_SLIP) :
		default: // defaults to BC_NO_SLIP
		{
			for(std::size_t i = istart; i <= iend; ++i) {
				dat->u_x[i][0] = -1.0*dat->u_x[i][1];
				dat->u_y[i][0] = 0;
			}
		} break;
	}
    
    // exchange first and last columns of u_x and u_y between all processes in order to treat internal cells bc correctly
    if(dat->rank > 0) {
        MPI_Send(&(dat->u_x[istart][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&(dat->u_x[istart - 1][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank - 1, 0, MPI_COMM_WORLD, &(dat->status));
        
        MPI_Send(&(dat->u_y[istart][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&(dat->u_y[istart - 1][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank - 1, 0, MPI_COMM_WORLD, &(dat->status));
    }
    
    if(dat->rank < dat->size - 1) {
        MPI_Send(&(dat->u_x[iend][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&(dat->u_x[iend + 1][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank + 1, 0, MPI_COMM_WORLD, &(dat->status));
        
        MPI_Send(&(dat->u_y[iend][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&(dat->u_y[iend + 1][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank + 1, 0, MPI_COMM_WORLD, &(dat->status));
    }

	// treat internal cell boundary conditions (no-slip only)
	for(std::size_t i = istart; i <= iend; i++) {
		for(std::size_t j = jstart; j <= jend; j++) {

			// The mask 0x000F filters the obstacle cells adjacent to  fluid cells
			if(dat->flag[i][j] & 0x000F) {
				switch(dat->flag[i][j]) {
					case B_N:
					{
						dat->u_y[i][j] = 0.0;
						dat->u_x[i][j] = -dat->u_x[i][j + 1];
						dat->u_x[i - 1][j] = -dat->u_x[i - 1][j + 1];
					} break;
					case B_E:
					{
						dat->u_x[i][j] = 0.0;
						dat->u_y[i][j] = -dat->u_y[i + 1][j];
						dat->u_y[i][j - 1] = -dat->u_y[i + 1][j - 1];
					} break;
					case B_S:
					{
						dat->u_y[i][j - 1] = 0.0;
						dat->u_x[i][j] = -dat->u_x[i][j - 1];
						dat->u_x[i - 1][j] = -dat->u_x[i - 1][j - 1];
					} break;
					case B_W:
					{
						dat->u_x[i - 1][j] = 0.0;
						dat->u_y[i][j] = -dat->u_y[i - 1][j];
						dat->u_y[i][j - 1] = -dat->u_y[i - 1][j - 1];
					} break;
					case B_NE:
					{
						dat->u_y[i][j] = 0.0;
						dat->u_x[i][j] = 0.0;
						dat->u_y[i][j - 1] = -dat->u_y[i + 1][j - 1];
						dat->u_x[i - 1][j] = -dat->u_x[i - 1][j + 1];
					} break;
					case B_SE:
					{
						dat->u_y[i][j - 1] = 0.0;
						dat->u_x[i][j] = 0.0;
						dat->u_y[i][j] = -dat->u_y[i + 1][j];
						dat->u_x[i - 1][j] = -dat->u_x[i - 1][j - 1];
					} break;
					case B_SW:
					{
						dat->u_y[i][j - 1] = 0.0;
						dat->u_x[i - 1][j] = 0.0;
						dat->u_y[i][j] = -dat->u_y[i - 1][j];
						dat->u_x[i][j] = -dat->u_x[i][j - 1];
					} break;
					case B_NW:
					{
						dat->u_y[i][j] = 0.0;
						dat->u_x[i - 1][j] = 0.0;
						dat->u_y[i][j - 1] = -dat->u_y[i - 1][j - 1];
						dat->u_x[i][j] = -dat->u_x[i][j + 1];
					} break;
					default: break;
				}
			}
		}
	}
}

void CFD_simulation::compute_F_G() {
    
    // aliases for easy access to boundary cells
	const unsigned& imax = dat->nmax[X];
	const unsigned& jmax = dat->nmax[Y];
    
    unsigned jstart = 1;
    unsigned jend = jmax;
    
    unsigned istart = dat->rank*imax/dat->size + 1;
    unsigned iend = (dat->rank + 1)*imax/dat->size;
    
    if(dat->rank < dat->size - 1) {
        for(std::size_t i = istart; i <= iend; ++i) {
            for(std::size_t j = jstart; j <= jend; ++j) {
                // only if both adjacent cells are fluid cells
                if((dat->flag[i][j] & C_F) && (dat->flag[i + 1][j] & C_F)) {
                    dat->F[i][j] = dat->u_x[i][j] + dat->delta_t * (
                        1.0 / dat->Re * (this->d2ux_dx2(i, j) + this->d2ux_dy2(i, j))
                        - this->dux2_dx(i, j) - this->duxuy_dy(i, j) + dat->g[X]);
                }
            }
        }
    }
    
    if(dat->rank == dat->size - 1) {
        for(std::size_t i = istart; i <= iend - 1; ++i) {
            for(std::size_t j = jstart; j <= jend; ++j) {
                // only if both adjacent cells are fluid cells
                if((dat->flag[i][j] & C_F) && (dat->flag[i + 1][j] & C_F)) {
                    dat->F[i][j] = dat->u_x[i][j] + dat->delta_t * (
                        1.0 / dat->Re * (this->d2ux_dx2(i, j) + this->d2ux_dy2(i, j))
                        - this->dux2_dx(i, j) - this->duxuy_dy(i, j) + dat->g[X]);
                }
            }
        }
    }
    
	for(std::size_t i = istart; i <= iend; ++i) {
		for(std::size_t j = jstart; j <= jend - 1; ++j) {
			// only if both adjacent cells are fluid cells
			if((dat->flag[i][j] & C_F) && (dat->flag[i][j + 1] & C_F)) {
				dat->G[i][j] = dat->u_y[i][j] + dat->delta_t * (
					1.0 / dat->Re * (this->d2uy_dx2(i, j) + this->d2uy_dy2(i, j))
					- this->duxuy_dx(i, j) - this->duy2_dy(i, j) + dat->g[Y]);
			}
		}
	}
}

void CFD_simulation::compute_Possion_RHS() {
	// aliases for easy access to boundary cells
	const unsigned& imax = dat->nmax[X];
	const unsigned& jmax = dat->nmax[Y];
    
    unsigned jstart = 1;
    unsigned jend = jmax;
    
    unsigned istart = dat->rank*imax/dat->size + 1;
    unsigned iend = (dat->rank + 1)*imax/dat->size;

	// set domain BC for F, G and p
    if(dat->rank == 0) {
        for(std::size_t j = jstart; j <= jend; ++j) {
            dat->p[0][j] = dat->p[1][j]; 				        // pressure west
            dat->F[0][j] = dat->u_x[0][j];				        // F west
        }
    }
    
    if(dat->rank == dat->size - 1) {
        for(std::size_t j = jstart; j <= jend; ++j) {
            dat->p[imax + 1][j] = dat->p[imax][j]; 		        // pressure east
            dat->F[imax][j] = dat->u_x[imax][j];		        // F east
        }
    }
    
	for(std::size_t i = istart; i <= iend; ++i) {
		dat->p[i][jmax + 1] = dat->p[i][jmax]; 		            // pressure north
        dat->G[i][jmax] = dat->u_y[i][jmax];		            // G north
    }
    
    for(std::size_t i = istart; i <= iend; ++i) {    
		dat->p[i][0] = dat->p[i][1]; 				            // pressure south
		dat->G[i][0] = dat->u_y[i][0];				            // G south
	}
    
    // send last column of F to next process in order to compute correct RHS
    if(dat->size > 1) {
        if(dat->rank == 0) {
            MPI_Send(&(dat->F[iend][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank + 1, 0, MPI_COMM_WORLD);
        }
        
        if(dat->rank > 0 && dat->rank < dat->size - 1) {
            MPI_Send(&(dat->F[iend][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&(dat->F[istart - 1][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank - 1, 0, MPI_COMM_WORLD, &(dat->status));
        }
        
        if(dat->rank == dat->size - 1) {
            MPI_Recv(&(dat->F[istart - 1][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank - 1, 0, MPI_COMM_WORLD, &(dat->status));
        }
    }

	// compute RHS
	for(std::size_t i = istart; i <= iend; ++i) {
		for(std::size_t j = jstart; j <= jend; ++j) {
			// only for fluid and non-surface cells
			if((dat->flag[i][j] & C_F) && (dat->flag[i][j] < 0x0100)) {
				dat->rhs[i][j] = 1.0 / dat->delta_t * ((dat->F[i][j] - dat->F[i - 1][j]) / dat->delta[X]
					+ (dat->G[i][j] - dat->G[i][j - 1]) / dat->delta[Y]);
			}
		}
	}
}

void CFD_simulation::solve_Poisson_Jacobi() {
	
    // local and global error
	double l_err, g_err;
    
    // iteration counter
    unsigned it;

	// alias for 1/delta_x^2 and 1/delta_y^2
	double idx2 = 1.0 / (dat->delta[X] * dat->delta[X]);
	double idy2 = 1.0 / (dat->delta[Y] * dat->delta[Y]);
    
    const unsigned& imax = dat->nmax[X];
    const unsigned& jmax = dat->nmax[Y];
    
    unsigned jstart = 1;
    unsigned jend = jmax;
    
    unsigned istart = dat->rank*imax/dat->size + 1;
    unsigned iend = (dat->rank + 1)*imax/dat->size;
    
	for(it = 1; it <= dat->itermax - 1; ++it) { 
        
		l_err = 0.0;
        g_err = 0.0;
        
        // exchange first and last column of new pressure between all domains in every iteration
        if(dat->rank > 0) {
            MPI_Send(&(dat->p[istart][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&(dat->p[istart - 1][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank - 1, 0, MPI_COMM_WORLD, &(dat->status));
        }
        
        if(dat->rank < dat->size - 1) {
            MPI_Send(&(dat->p[iend][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&(dat->p[iend + 1][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank + 1, 0, MPI_COMM_WORLD, &(dat->status));
        }
        
		for(std::size_t i = istart; i <= iend; ++i) {
			for(std::size_t j = jstart; j <= jend; ++j) {
				
                // pure fluid cells in the inside with no connections to boundaries
				if(dat->flag[i][j] == 0x001F) {
					dat->p_new[i][j] = (1.0 - dat->omega) * dat->p[i][j] + 
                        dat->omega / (2.0*idx2 + 2.0*idy2) * ((dat->p[i + 1][j] +
                        dat->p[i - 1][j])*idx2 + (dat->p[i][j + 1] + 
                        dat->p[i][j - 1])*idy2 - dat->rhs[i][j]);
					// determine local error by using the maximum norm
					if(fabs(dat->p[i][j] - dat->p_new[i][j]) > l_err) {
						l_err = fabs(dat->p[i][j] - dat->p_new[i][j]);
					}
				} 
                
                // compute fluid cells with connections to some boundaries
				if((dat->flag[i][j] & C_F) && (dat->flag[i][j] < 0x0100)) {
					dat->p_new[i][j] = (1.0 - dat->omega) * dat->p[i][j] + 
                        dat->omega / ((eps_E + eps_W)*idx2 + (eps_N + eps_S)*idy2) *
						((eps_E*dat->p[i + 1][j] + eps_W*dat->p[i - 1][j])*idx2 + 
                        (eps_N*dat->p[i][j + 1] + eps_S*dat->p[i][j - 1])*idy2 -
                        dat->rhs[i][j]);
					// determine local error by using the maximum norm
					if(fabs(dat->p[i][j] - dat->p_new[i][j]) > l_err) {
						l_err = fabs(dat->p[i][j] - dat->p_new[i][j]);
					}
				}
                
			}
		}
        
        // swap the interior points
        for(std::size_t i = istart; i <= iend; ++i) {
			for(std::size_t j = jstart; j <= jend; ++j) {
                dat->p[i][j] = dat->p_new[i][j];
            }
        }
        
        // wait for all processes to compute local error in their domain
        MPI_Barrier(MPI_COMM_WORLD);
        
        // assign the highest local error to global error
        MPI_Allreduce(&l_err, &g_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
             
        // terminate iteration if global error satisfies tolerance
		if(g_err < dat->eps) {
            break;
        }
	}
    // debug info
    if(dat->rank == 0) {
        std::cout << "solved Poisson eq. for t = " << dat->t << " using " << it << " iterations (max err = " << g_err << ")" << std::endl;
    }
}

void CFD_simulation::compute_u_next() {
       
    const unsigned& imax = dat->nmax[X];
    const unsigned& jmax = dat->nmax[Y];
    
    unsigned jstart = 1;
    unsigned jend = jmax;
    
    unsigned istart = dat->rank*imax/dat->size + 1;
    unsigned iend = (dat->rank + 1)*imax/dat->size;
    
	std::array<double, 2> u_max_abs_loc = {0.0, 0.0};
    
    // send first column of new p to the previous process in order to compute correct u next
    if(dat->size > 1) {
        if(dat->rank == 0) {
            MPI_Recv(&(dat->p[iend + 1][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank + 1, 0, MPI_COMM_WORLD, &(dat->status));
        }
        
        if(dat->rank > 0 && dat->rank < dat->size - 1) {
            MPI_Send(&(dat->p[istart][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&(dat->p[iend + 1][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank + 1, 0, MPI_COMM_WORLD, &(dat->status));
        }
        
        if(dat->rank == dat->size - 1) {
            MPI_Send(&(dat->p[istart][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank - 1, 0, MPI_COMM_WORLD);
        }
    }
    
    if(dat->rank < dat->size - 1) {
        for(std::size_t i = istart; i <= iend; ++i) {
            for(std::size_t j = jstart; j <= jend; ++j) {
                // only if both adjacent cells are fluid cells
                if((dat->flag[i][j] & C_F) && (dat->flag[i + 1][j] & C_F)) {
                    dat->u_x[i][j] = dat->F[i][j] - dat->delta_t / dat->delta[X] * (dat->p[i + 1][j] - dat->p[i][j]);
                    if(fabs(dat->u_x[i][j]) > u_max_abs_loc[X]) {
                        u_max_abs_loc[X] = fabs(dat->u_x[i][j]);
                    }
                }
            }
        }
    }
    
    // last process cannot acces flag[iend+1][j]
    if(dat->rank == dat->size - 1) {
        for(std::size_t i = istart; i <= iend - 1; ++i) {
            for(std::size_t j = jstart; j <= jend; ++j) {
                // only if both adjacent cells are fluid cells
                if((dat->flag[i][j] & C_F) && (dat->flag[i + 1][j] & C_F)) {
                    dat->u_x[i][j] = dat->F[i][j] - dat->delta_t / dat->delta[X] * (dat->p[i + 1][j] - dat->p[i][j]);
                    if(fabs(dat->u_x[i][j]) > u_max_abs_loc[X]) {
                        u_max_abs_loc[X] = fabs(dat->u_x[i][j]);
                    }
                }
            }
        }         
    }

	for(std::size_t i = istart; i <= iend; ++i) {
		for(std::size_t j = jstart; j <= jend - 1; ++j) {
			// only if both adjacent cells are fluid cells
			if((dat->flag[i][j] & C_F) && (dat->flag[i][j + 1] & C_F)) {
				dat->u_y[i][j] = dat->G[i][j] - dat->delta_t / dat->delta[Y] * (dat->p[i][j + 1] - dat->p[i][j]);
				if(fabs(dat->u_y[i][j]) > u_max_abs_loc[Y]) {
                    u_max_abs_loc[Y] = fabs(dat->u_y[i][j]);
                }
			}
		}
	}
    
    // wait for all processes to compute maximum absolute velocities in their domain
    MPI_Barrier(MPI_COMM_WORLD);
    
    // get the maximum absolute velocities from all processes to compute correct gamma value for all processes  
    MPI_Allreduce(&(u_max_abs_loc[X]), &(dat->u_max_abs[X]), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&(u_max_abs_loc[Y]), &(dat->u_max_abs[Y]), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    // send last u_x column to next process for visualization purposes 
    if(dat->size > 1) {
        if(dat->rank == 0) {
            MPI_Send(&(dat->u_x[iend][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank + 1, 0, MPI_COMM_WORLD);
        }
        
        if(dat->rank > 0 && dat->rank < dat->size - 1) {
            MPI_Send(&(dat->u_x[iend][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&(dat->u_x[istart - 1][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank - 1, 0, MPI_COMM_WORLD, &(dat->status));
        }
        
        if(dat->rank == dat->size - 1) {
            MPI_Recv(&(dat->u_x[istart - 1][0]), dat->nmax[Y] + 2, MPI_DOUBLE, dat->rank - 1, 0, MPI_COMM_WORLD, &(dat->status));
        }
    }
}

void CFD_simulation::compute_step() {
    this->adapt_delta_t();
    this->update_domain();
    this->compute_F_G();
    this->compute_Possion_RHS();
    this->solve_Poisson_Jacobi();
    this->compute_u_next();
    // advance t 
    dat->t += dat->delta_t;
}
