#include "inc/simulation_data.hpp"

simulation_data::simulation_data(const std::array<double, 2> _length, const std::array<unsigned int, 2> _nmax, int _rank, int _size)
	:length(_length), nmax(_nmax), rank(_rank), size(_size) {
	
	delta = {
		length[X] / (double)nmax[X],
		length[Y] / (double)nmax[Y]
	};
	
	// default values
	t = 0.0;
	t_end = 10.0;
	delta_t = 1e-2;
	delta_t_is_maximal_value = false;
	tau = 0.75;
	itermax = 10000;
	eps = 1e-1;
	omega = 1.0;
	gamma = 0.0;
	Re = 100;
	g = {0.0, 0.0};
	u_inflow = {0.0, 0.0, 0.0, 0.0};
	u_max_abs = {0.0, 0.0};
	
	domain_boundary[D_EAST] = BC_NO_SLIP;
	domain_boundary[D_WEST] = BC_NO_SLIP;
	domain_boundary[D_NORTH] = BC_NO_SLIP;
	domain_boundary[D_SOUTH] = BC_NO_SLIP;
	
	u_x.resize(nmax[X] + 2);
	u_y.resize(nmax[X] + 2);
	p.resize(nmax[X] + 2);
    p_new.resize(nmax[X] + 2);
	rhs.resize(nmax[X] + 2);
	F.resize(nmax[X] + 2);
	G.resize(nmax[X] + 2);
	flag.resize(nmax[X] + 2);
	
	for(std::size_t i = 0; i < nmax[X] + 2; ++i) {
        
		u_x[i].resize(nmax[Y] + 2);
		u_y[i].resize(nmax[Y] + 2);
		p[i].resize(nmax[Y] + 2);
        p_new[i].resize(nmax[Y] + 2);
		rhs[i].resize(nmax[Y] + 2);
		F[i].resize(nmax[Y] + 2);
		G[i].resize(nmax[Y] + 2);
		flag[i].resize(nmax[Y] + 2);
        
		for(std::size_t j = 0; j < nmax[Y] + 2; ++j) {
			flag[i][j] = C_F;
		}
	}
	
	return;
}

void simulation_data::initial_conditions(const double& u_x_init, const double& u_y_init, const double& p_init) {
	u_max_abs = {fabs(u_x_init), fabs(u_y_init)};
	for(std::size_t i = 0; i < nmax[X] + 2; ++i) {
		for(std::size_t j = 0; j < nmax[Y] + 2; ++j) {
			u_x[i][j] = u_x_init;
			u_y[i][j] = u_y_init;
			p[i][j] = p_init;
		}
	}
	return;
}

bool simulation_data::boundary_conditions(std::vector<std::vector<int>> obstacle_matrix) {

	// mark all the boundaries in the flag field
	for(std::size_t i = 0; i <= this->nmax[X] + 1; ++i) {
		flag[i][0] = C_B;
		flag[i][this->nmax[Y] + 1] = C_B;
	}
	for(std::size_t j = 1; j <= this->nmax[Y]; ++j) {
		flag[0][j] = C_B;
		flag[this->nmax[X] + 1][j] = C_B;
	}

	int entry = 0;
	for(std::size_t i = 1; i <= this->nmax[X]; ++i) {
		for(std::size_t j = 1; j <= this->nmax[Y]; ++j) {
			entry = obstacle_matrix[i - 1][j - 1];
			// only values of C_B and C_F are accepted
			// if other values are present, they will be set to C_B by default
			// a warning will be issued
			if(!(entry == C_F || entry == C_B)) {
                if(rank == 0) {
                    std::cout << "warning: entry " << i - 1 << "," << j - 1 << " has a value of "
                        << entry << " which is non C_B or C_F value. C_B is assumed." << std::endl;
                }
				entry = C_B;
			}
			flag[i][j] = entry;
		}
	}
	
	// adjust flag field to reflect different fluid cell states with boundary cells
	int count_boundary_cells = 0;
	for(std::size_t i = 1; i <= this->nmax[X]; ++i) {
		for(std::size_t j = 1; j <= this->nmax[Y]; ++j) {
			// check for boundary cells
			if(!(flag[i][j] & C_F)) count_boundary_cells++;
			// adjust flag field
			flag[i][j] += ((flag[i - 1][j] & C_F)*B_W + (flag[i + 1][j] & C_F)*B_E +
						  (flag[i][j - 1] & C_F)*B_S + (flag[i][j + 1] & C_F)*B_N) / C_F;
			// check if cell got an invalid flag of an invalid obstacle cell
			switch(flag[i][j]) {
				case 0x0003:
				case 0x0007:
				case 0x000b:
				case 0x000c:
				case 0x000d:
				case 0x000e:
				case 0x000f:
				{
                    if(rank == 0) {
                        std::cout << "warning: illegal obstacle cell [" << i << "][" << j << "] detected." << std::endl;
                    }
                    return false;
				} break;
			}
		}
	}
	// debug info
	if(count_boundary_cells) {
        if(rank == 0) {
            std::cout << "number of cells: " << nmax[X] * nmax[Y] << " boundary cells: " << count_boundary_cells << std::endl;
        }
    }
	return true;
}

void simulation_data::export_vtk(const std::string filename) const {
    
    const unsigned& imax = this->nmax[X];
    const unsigned& jmax = this->nmax[Y];
    
    unsigned jstart = 1;
    unsigned jend = jmax;
    
    unsigned istart = this->rank*imax/this->size + 1;
    unsigned iend = (this->rank + 1)*imax/this->size;

	std::ofstream output_file(filename, std::ios::out | std::ios::app);
	if(!output_file) {
		std::cerr << "error: couldn't write file " << filename << std::endl;
		return;
	}
    
    if(this->rank == 0) {
        output_file
            << "# vtk DataFile Version 2.0" << std::endl
            << "MESH" << std::endl
            << "ASCII" << std::endl
            << "DATASET RECTILINEAR_GRID" << std::endl
            << "DIMENSIONS " << this->nmax[X] + 1 << " " << this->nmax[Y] + 1 << " " << 1 << std::endl;

        output_file
            << "X_COORDINATES " << this->nmax[X] + 1 << " float" << std::endl;
        double coord = 0.0;
        for(std::size_t i = 0; i <= this->nmax[X]; ++i) {
            output_file << coord << std::endl;
            coord += this->delta[X];
        }

        output_file
            << "Y_COORDINATES " << this->nmax[Y] + 1 << " float" << std::endl;
        coord = 0.0;
        for(std::size_t j = 0; j <= this->nmax[Y]; ++j) {
            output_file << coord << "\n";
            coord += this->delta[Y];
        }

        output_file
            << "Z_COORDINATES " << 1 << " float" << std::endl
            << 0.0 << std::endl;

        output_file
            << "CELL_DATA " << this->nmax[X] * this->nmax[Y] << std::endl
            << "VECTORS velocities float" << std::endl;
    }
    
    int count1;
	for(std::size_t j = jstart; j <= jend; ++j) {
        count1 = 0;
        while(count1 < this->size) {
            if(this->rank == count1) {
                for(std::size_t i = istart; i <= iend; ++i) {
                        output_file
                            << 0.5*(this->u_x[i - 1][j] + this->u_x[i][j]) << " "
                            << 0.5*(this->u_y[i][j - 1] + this->u_y[i][j]) << " "
                            << 0.0 << std::endl;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            count1++;
        }
	}
    
    if(this->rank == 0) {
        output_file
            << "SCALARS pressure float" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
    }
    
    int count2;
	for(std::size_t j = jstart; j <= jend; ++j) {
        count2 = 0;
        while(count2 < this->size) {
            if(this->rank == count2) {
                for(std::size_t i = istart; i <= iend; ++i) {
                    output_file << this->p[i][j] << std::endl;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            count2++;
		}
	}
    
    if(this->rank == 0) {
        output_file
            << "SCALARS flag int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
        for(std::size_t j = 1; j <= this->nmax[Y]; ++j) {
            for(std::size_t i = 1; i <= this->nmax[X]; ++i) {
                output_file << this->flag[i][j] << std::endl;
            }
        }
    }

	output_file.close();
	return;
}