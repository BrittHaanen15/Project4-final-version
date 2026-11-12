//Unit test to check that the init and padding of my class works as expected

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <chrono>
#include <algorithm>
#include "utils.hpp"
#include "eigen-5.0.1/eigen-5.0.1/Eigen/Dense"
#include "eigen-5.0.1/eigen-5.0.1/Eigen/Eigenvalues"

#include "Ising2D.hpp"

int main()
{
    int lattice_size = 2;
    double temp = 1.0;

    //For this unit test I hard coded the values to put in the lattice.
    // //Fully control how it's set for now
    // m1(0,0) = 1;
    // m1(0,1) = 1;
    // m1(1,0) = 1;
    // m1(1,1) = 1;
    //Also tried:
    // m1(0,0) = -1;
    // m1(0,1) = 1;
    // m1(1,0) = -1;
    // m1(1,1) = 1;
    //////////////////////////////////////////////

    //Initialise
    Ising2D model = Ising2D(lattice_size, temp);

    model.print();

    //Calculate and print values
    //Expect (-16)/4=-4 for epsilon and (4)/4=1 for m (first case)
    //Expect 0 and 0 for second case

    std::cout << "epsilon: " << model.calculate_epsilon() << std::endl;
    std::cout << "|m|: " << model.calculate_mag() << std::endl;

    return 0;
}
