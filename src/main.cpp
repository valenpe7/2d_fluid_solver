#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include <mpi.h>

#include "inc/simulation_data.hpp"
#include "inc/cfd_simulation.hpp"
#include "inc/lib/lodepng/lodepng.h"

typedef std::chrono::high_resolution_clock Clock;

void run_simulation(const std::string filename, const std::string output_path, int np, int id) {

    unsigned size_x, size_y;
    std::vector<unsigned char> image;
    std::vector<std::vector<int>> obstacle;
    
    std::array<double, 2> length;
    std::array<unsigned, 2> nmax;
    
	// master reads png image
    if (id == 0) {
        unsigned error = lodepng::decode(image, size_x, size_y, filename);
        if(error) {
            std::cerr << "decoder error" << error << ": " << lodepng_error_text(error) << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
            return;
        }
        if (size_x % np != 0) {
            std::cout << "error: number of cells has to be a multiple of number of processors " << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
            return;
        }
    }
    
    // broadcast size of image to all processes
    MPI_Bcast(&size_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&size_y, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // resize obstacle matrix
    obstacle.resize(size_x);
	for(std::size_t i = 0; i < size_x; ++i) {
		obstacle[i].resize(size_y);
    }
    
    // master computes obstacle matrix from image
    if (id == 0) {
        int entry = 0;
        for(std::size_t i = 0; i < size_x; ++i) {
            for(std::size_t j = 0; j < size_y; ++j) {
                entry = image[4 * ((size_y - j - 1) * size_x + i) + 3];
                if(!entry) {
                    obstacle[i][j] = C_F;
                } else {
                    obstacle[i][j] = C_B;
                }
            }
        }
    }
    
    // broadcast obstacle matrix to all processes
    for(std::size_t i = 0; i < size_x; ++i) {
        MPI_Bcast(&(obstacle[i][0]), size_y, MPI_INT, 0, MPI_COMM_WORLD);
    }
    
    // set length from image and scale it by 40 (good value for the channel)
    length = {size_x / 40.0, size_y / 40.0};

    // set number of cells according to image size
    nmax = {size_x, size_y};

	// create data structure necessary for simulation
	simulation_data sim_data(length, nmax, id, np);

	// initialize velocity and pressure fields with initial values for channel flow
	double u_x_init = 0.0, u_y_init = 0.0, p_init = 0.0;
	sim_data.initial_conditions(u_x_init, u_y_init, p_init);
    
    // set interior obstacles
	sim_data.boundary_conditions(obstacle);

	// set CFD channel properties
	sim_data.Re = 100.0;
	sim_data.itermax = 10000;
	sim_data.t_end = 20.0;
	sim_data.delta_t = 1e-3;
	sim_data.eps = 1e-1;
	sim_data.omega = 0.67; // optimal value

    // inflow in the channel (define at the west boundary)
	sim_data.domain_boundary[D_NORTH] = BC_INFLOW;
    sim_data.u_inflow[D_NORTH] = 1.5;
    
    // outflow of the channel (define at the south boundary)
	sim_data.domain_boundary[D_SOUTH] = BC_OUTFLOW;
    
    // sides of the channel
	sim_data.domain_boundary[D_WEST] = BC_NO_SLIP;
	sim_data.domain_boundary[D_EAST] = BC_NO_SLIP;

	// create the simulation class and link the data
	CFD_simulation simulation(&sim_data);

	// visualization interval
	double t_vis = 0.5;
	double t_next_vis = 0.0;
    
    // file counter
	int count = 0;
    std::string output;

	// main CFD loop
	while(sim_data.t < sim_data.t_end) {
        t_next_vis += t_vis;
		while(sim_data.t < t_next_vis) {
			simulation.compute_step();
		}
		// export data
        std::ostringstream stream;
        stream << output_path << "/vis_" << std::setw(4) << std::setfill('0') << count++ << ".vtk";
        output = stream.str();
        // debug info:
        // if(sim_data.rank == 0) {
        //      std::cout << "writing data into file: " << output << std::endl;
        // }
        // for performance tests without exporting the data!!
        // sim_data.export_vtk(output);
        // MPI_Barrier(MPI_COMM_WORLD);
	}

	return;
}

int main(int argc, char** argv) {
	if (argc != 3) {
		std::cout << "usage: ./application_name.exe input_file.png output_folder" << std::endl;
        return -1;
    }
    
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// debug: performance test
	auto t1 = Clock::now();
    
	run_simulation(argv[1], argv[2], size, rank);
    
	//debug: performance test
	auto t2 = Clock::now();
    
    MPI_Finalize();
    
    //debug: performance test
    if(rank == 0) {
    	std::cout << "number of cores: " << size << std::endl;
    	std::cout << "duration: " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " s" << std::endl;
    }

    return 0;
}
