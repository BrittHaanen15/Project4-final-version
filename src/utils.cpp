#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <chrono>
#include <random>
#include "utils.hpp"
#include "../include/eigen-5.0.1/eigen-5.0.1/Eigen/Dense"
#include "../include/eigen-5.0.1/eigen-5.0.1/Eigen/Eigenvalues"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Make a vector according to start, stop, step
std::vector<double> make_a_range(double start, double stop, double step, bool verbose)
{
    //Calc number of steps based on start, stop and step size
    int stepsnumber = static_cast<int>((stop-start)/step);

    //Preall vector of this size
    std::vector<double> output(stepsnumber);

    for (int i = 0; i < stepsnumber; i++)
    {
        output[i] = start + i*step;
    }

    //Reality check: print generated vector if needed
    if (verbose == true)
    {
        for (const auto& element : output) 
        {
            std::cout << element << " ";
        }
        std::cout << "\n";
    }
    
   return output;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



