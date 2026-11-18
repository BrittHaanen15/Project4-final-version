//Search for phase transitions, over different lattice sizes. Follows the parallellisation structure of Skeleton_parallellised.cpp by threading over temperatures
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <random>
#include "utils.hpp"
#include <ctime>
#include <omp.h>
#include "eigen-5.0.1/eigen-5.0.1/Eigen/Dense"
#include "eigen-5.0.1/eigen-5.0.1/Eigen/Eigenvalues"

#include "Ising2D.hpp"

int main()
{
    //Lattice sizes to try
    std::vector<double> lattices = make_a_range(40, 120, 20, false);

    for (int j = 0; j<static_cast<int>(lattices.size()); j++) //Run 5 timing tests
    {
        //Init variables
        int lattice_size = lattices[j];
        int n_cycles = 1000000;
        int steps_in_cycle = 10;
        bool ordered_state = false;

        //Temperatures to try
        std::vector<double> temperatures = make_a_range(2.23, 2.25, 0.001, true);

        #ifdef _OPENMP
        {
            //Parallellise over several temperatures (level 1 as given at https://anderkve.github.io/FYS3150/book/projects/project4.html#problem-7)
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(temperatures.size()); i++)
            {
                double temp = temperatures[i];

                //File for writing results
                std::string filename = "Outputs/VERYFINESEARCH/PhaseTrans_L" + std::to_string(lattice_size) + "_T" + std::to_string(temp) + "_unordered_n10e6.txt";
                
                std::ofstream outFile(filename);
                outFile << "epsilon, abs(m)\n";

                // Get a seed from the system clock
                //Adapted from code by Anders Kvellestad at https://github.com/anderkve/FYS3150/tree/master/code_examples/random_number_generation
                unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();

                //Initialise
                Ising2D model = Ising2D(lattice_size, temp, seed, ordered_state);

                for (int i=0; i<n_cycles; i++)
                {
                    for (int j=0; j<steps_in_cycle; j++)
                    {
                        model.flip_a_spin();
                    }
                    
                    //Calculate
                    double eps = model.calculate_epsilon();
                    double mag = model.calculate_mag();

                    //Write
                    outFile << eps << ","
                            << mag << "\n";

                }
            }
        }
        #else
        {
            for (int i = 0; i < static_cast<int>(temperatures.size()); i++)
                {
                    double temp = temperatures[i];

                    //File for writing results
                    std::string filename = "Outputs/VERYFINESEARCH/PhaseTrans_L" + std::to_string(lattice_size) + "_T" + std::to_string(temp) + "_unordered_n10e6.txt";
                    
                    std::ofstream outFile(filename);
                    outFile << "epsilon, abs(m)\n";

                    // Get a seed from the system clock
                    //Adapted from code by Anders Kvellestad at https://github.com/anderkve/FYS3150/tree/master/code_examples/random_number_generation
                    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();

                    //Initialise
                    Ising2D model = Ising2D(lattice_size, temp, seed, ordered_state);

                    for (int i=0; i<n_cycles; i++)
                    {
                        for (int j=0; j<steps_in_cycle; j++)
                        {
                            model.flip_a_spin();
                        }
                        
                        //Calculatecd
                        double eps = model.calculate_epsilon();
                        double mag = model.calculate_mag();

                        //Write
                        outFile << eps << ","
                                << mag << "\n";

                    }
                }
        }
        #endif

        std::cout << "Finished lattice size " << lattice_size << std::endl;
    }

    return 0;
}

