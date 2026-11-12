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
    //File for writing results
    std::ofstream outFile("Outputs/Timingtest_YES.txt");
    outFile << "time(s)\n";
    //Start timing
    std::clock_t t1 = std::clock();

    for (int j = 0; j<5; j++) //Run 10 timing tests
    {
        

        //Init variables
        int lattice_size = 20;
        int n_cycles = 1000000;
        int steps_in_cycle = 10;
        bool ordered_state = true;

        //Temperatures to try
        std::vector<double> temperatures = make_a_range(1.0, 4.2, 0.2, false);

        #ifdef _OPENMP
        {
            //Parallellise over several temperatures (level 1 as given at https://anderkve.github.io/FYS3150/book/projects/project4.html#problem-7)
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(temperatures.size()); i++)
            {
                double temp = temperatures[i];

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

                }
            }
        }
        #else
        {
            for (int i = 0; i < static_cast<int>(temperatures.size()); i++)
                {
                    double temp = temperatures[i];

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

                    }
                }
        }
        #endif
    }

    // Stop measuring time
    std::clock_t t2 = std::clock();

    // Calculate the elapsed time.
    double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;

    //Write
    outFile << duration_seconds << "\n";

    return 0;
}
