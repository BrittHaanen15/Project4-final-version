//2D Ising model class header file

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "eigen-5.0.1/eigen-5.0.1/Eigen/Dense"
#include "eigen-5.0.1/eigen-5.0.1/Eigen/Eigenvalues"
#ifndef __Ising2D_hpp__
#define __Ising2D_hpp__

class Ising2D
{
private:
  // Nothing private here so far, right?

public:
  //Member variables////////////////////
  //Lattice size L, temperature T, given at init
  int L_;
  double T_;
  //Lattice and padded lattice
  Eigen::MatrixXd lattice_, P_lattice_;
  //Random number generator
  std::uniform_int_distribution<int> my_spin_pdf_;
  std::uniform_int_distribution<int> my_0L_pdf_;
  std::uniform_real_distribution<double> my_01_pdf_;
  std::mt19937 generator_;
  //Constants
  double beta;

  // Constructor
  Ising2D(int Lattice_size, double Temperature, int seed = 2026);

  //Initialise lattice of specified size L_
  Eigen::MatrixXd set_lattice(int size);

  //Pad lattice for calculations using periodic boundary conditions
  Eigen::MatrixXd pad_lattice(Eigen::MatrixXd lat);

  //Set up random nr gen
  std::mt19937 set_up_generator(int seed);

  //Print lattices for checking
  void print();

  //Calculate energy per spin
  double calculate_epsilon();

  //Calculate magnetisation per spin;
  double calculate_mag();

  //Flip a spin randomly
  void flip_a_spin();

  //Calculate energy of one spin 
  double one_spin_energy(int row, int col, double spin_value);

  //return boltzmann factor for given energy change
  double return_boltzmann_factor(double energy_change);
};

#endif
