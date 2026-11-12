//2D ising model source file
#include "Ising2D.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include "eigen-5.0.1/eigen-5.0.1/Eigen/Dense"
#include "eigen-5.0.1/eigen-5.0.1/Eigen/Eigenvalues"

Ising2D::Ising2D(int Lattice_size, double Temperature, int seed)
{
    //Assign member variables from the given inputs
    L_ = Lattice_size;
    T_ = Temperature;

    //Random nr distributions
    //call as my_01_pdf_(generator_)
    //This one is for 0/1 (so for setting spins)
    my_spin_pdf_ = std::uniform_int_distribution<int>(0, 1);
    //This one is for selecting a spin to flip so for proposing a new candidate state
    my_0L_pdf_ = std::uniform_int_distribution<int>(0, (L_-1)); 
    //This one is for accepting or rejecting, generate a double on 0-1
    my_01_pdf_ = std::uniform_real_distribution<double>(0,1);

    //Construct generator
    generator_ = Ising2D::set_up_generator(seed);

    //Initialise lattice according to specified size
    lattice_ = Ising2D::set_lattice(L_);

    //Pad lattice for calculations
    P_lattice_ = Ising2D::pad_lattice(lattice_);

    //Calculate inverse temp
    beta = 1/T_;
}



//Initialise lattice of specified size L_
Eigen::MatrixXd Ising2D::set_lattice(int size)
{
    Eigen::MatrixXd m1 = Eigen::MatrixXd::Zero(size, size);

    // Hard-coded like this for the unit testing///
    // //Fully control how it's set for now
    // m1(0,0) = 0;
    // m1(0,1) = 1;
    // m1(1,0) = 2;
    // m1(1,1) = 3;
    //////////////////////////////////////////////

    //For normal algorithm: go with random initialisation
    for(int i=0; i < L_; i++)
    {
        for(int j=0; j < L_; j++)
        {
            int spin = my_spin_pdf_(generator_);
            if (spin==0)
            {
                m1(i,j) = -1;
            }
            else
            {
                m1(i,j) = spin;
            }
            
        }
    }

    return m1;
}



//Pad lattice for calculations using periodic boundary conditions
Eigen::MatrixXd Ising2D::pad_lattice(Eigen::MatrixXd lat)
{
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero((lat.rows()+2), (lat.cols() +2));
    out.block(1, 1, lat.rows(), lat.cols()) = lat;

    //row 0 of padded becomes last row of lat
    out.row(0).segment(1, lat.cols()) = lat.row(L_ - 1); 
    //last row of padded becomes first row of lat
    out.row(out.rows() - 1).segment(1, lat.cols()) = lat.row(0); 
    //col 0 of padded becomes last col of lat
    out.col(0).segment(1, lat.rows()) = lat.col(L_ - 1);   
    //last col of padded gets first col of lat
    out.col(out.cols() - 1).segment(1, lat.rows()) = lat.col(0);

    return out;
}



//Set up random nr gen
std::mt19937 Ising2D::set_up_generator(int seed)
{
  //Adapted from code by Anders Kvellestad at https://github.com/anderkve/FYS3150/tree/master/code_examples/random_number_generation
  // Construct a Mersenne Twister 19937 random number generator
  std::mt19937 generator;

  // Seed it with our seed
  generator.seed(seed);

  return generator;
}



//Print lattices for checking
void Ising2D::print()
{
    std::cout << lattice_ << std::endl;

    std::cout << "####################" << std::endl;

    std::cout << P_lattice_ << std::endl;
}



//Calculate energy per spin
double Ising2D::calculate_epsilon()
{
    double epsilon = 0;

    //Loop through elements of padded lattice
    //Only use inner elements
    for(int i=1; i < L_+1; i++)
    {
        for(int j=1; j < L_+1; j++)
        {
            epsilon += -P_lattice_(i,j)*P_lattice_(i+1,j);
            epsilon += -P_lattice_(i,j)*P_lattice_(i,j+1);
        }
    }


    //Calculate per spin
    epsilon = epsilon/static_cast<double>(lattice_.size());

    return epsilon;
}



//Calculate magnetisation per spin;
double Ising2D::calculate_mag()
{
    double mag = 0;
    //Loop through elements of padded lattice
    //Only use inner elements
    for(int i=1; i < L_+1; i++)
    {
        for(int j=1; j < L_+1; j++)
        {
            mag += -P_lattice_(i,j);
        }
    }

    //Calculate per spin
    mag = std::abs(mag)/static_cast<double>(lattice_.size());

    return mag;
}



// Randomly flip a spin
void Ising2D::flip_a_spin()
{
    //Pick an index
    int spin_row = my_0L_pdf_(generator_);
    int spin_col = my_0L_pdf_(generator_);

    // std::cout << "Indexes selected: " << spin_row << "   " << spin_col << std::endl;

    //Convert to padded lattice indexes
    spin_row += 1;
    spin_col += 1;
    //retrieve value
    int spin = P_lattice_(spin_row, spin_col);

    //Compute ratio
    //Can simplify the fraction (p(B))/(p(A)) by taking out the 1/Z and using the rule that fraction of exp = keep base with exp upper-lower
    //So the ratio becomes e^(-beta*deltaE)
    //deltaE can only take 5 possible values like we talked about for problem 2
    double E_after = Ising2D::one_spin_energy(spin_row, spin_col, (-1)*spin);
    double E_before = Ising2D::one_spin_energy(spin_row, spin_col, spin);
    double delta_E = E_after - E_before;
    // std::cout << "Energy change: " << delta_E << std::endl;
    double ratio = Ising2D::return_boltzmann_factor(delta_E);

    // std::cout << "ratio " << ratio << std::endl;

    //Accept or reject
    double r = my_01_pdf_(generator_);
    // std::cout << "random " << r << std::endl;
    if (r<=ratio)
    {
        //Accept
        // std::cout << "Accepted" << std::endl;

        //Flip the spin
        lattice_(spin_row-1, spin_col-1) = (-1)*lattice_(spin_row-1, spin_col-1);
        P_lattice_ = Ising2D::pad_lattice(lattice_);
    }
    else
    {
        //Reject
        // std::cout << "Rejected" << std::endl;

        //Spin will not be flipped
    }

}



//Calculate energy of one spin
double Ising2D::one_spin_energy(int row, int col, double spin_value)
{
    double e = 0;

    //Circle around indicated spin to count up the numbers  
    e += -spin_value*P_lattice_(row-1,col);
    e += -spin_value*P_lattice_(row+1,col);
    e += -spin_value*P_lattice_(row,col-1);
    e += -spin_value*P_lattice_(row,col+1);


    return e;
}


double Ising2D::return_boltzmann_factor(double energy_change)
{
    return std::exp(-beta*energy_change);
}