#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <random>
#include "utils.hpp"
#include "eigen-5.0.1/eigen-5.0.1/Eigen/Dense"
#include "eigen-5.0.1/eigen-5.0.1/Eigen/Eigenvalues"

#include "Ising2D.hpp"

int main()
{
    int lattice_size = 20;
    double temp = 2.4;
    int n_cycles = 1000000;
    int steps_in_cycle = 10;
    bool ordered_state = true;
    //File for writing results
    std::string filename = "Outputs/L20_T" + std::to_string(temp) + "_ordered_n10e6.txt";
    
    std::ofstream outFile(filename);
    outFile << "epsilon, abs(m)\n";

    // Get a seed from the system clock
    //Adapted from code by Anders Kvellestad at https://github.com/anderkve/FYS3150/tree/master/code_examples/random_number_generation
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();

    //Initialise
    Ising2D model = Ising2D(lattice_size, temp, seed, ordered_state);
    model.print();

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

    return 0;
}
