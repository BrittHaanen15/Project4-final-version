//Unit test to check that the randomisation and flipping of spins works as expected

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
    int lattice_size = 2;
    double temp = 1.0;
    int n_cycles = 5;
    int steps_in_cycle = 1;
    int seed = 2007; //Played around with this a bit until I found one that actually causes changes
    //File for writing results
    std::ofstream outFile("Outputs/UnitTest_Randomisation.txt");
    outFile << "epsilon, abs(m)\n";

    //Initialise
    Ising2D model = Ising2D(lattice_size, temp, seed);
    model.print();

    for (int i=0; i<n_cycles; i++)
    {
        for (int j=0; j<steps_in_cycle; j++)
        {
            model.flip_a_spin();
            model.print();
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
