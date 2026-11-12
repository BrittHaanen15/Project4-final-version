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
    //See also Ising2D.cpp
    // Hard-coded like this for the unit testing///
    // //Fully control how it's set for now
    // m1(0,0) = 0;
    // m1(0,1) = 1;
    // m1(1,0) = 2;
    // m1(1,1) = 3;
    //////////////////////////////////////////////

    //Initialise
    Ising2D model = Ising2D(lattice_size, temp);

    //Print
    model.print();

    return 0;
}
