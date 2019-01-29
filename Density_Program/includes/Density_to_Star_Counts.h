#ifndef __DENSITY_TO_STAR_COUNTS_H__
#define __DENSITY_TO_STAR_COUNTS_H__
#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <stdlib.h>
#include <random>
#include <time.h>
#include <cstdlib>
#include <cassert>

#include "synchronous_algorithms/parameter_sweep.hxx"
#include "synchronous_algorithms/synchronous_gradient_descent.hxx"
#include "synchronous_algorithms/synchronous_newton_method.hxx"

//from undvc_common
#include "arguments.hxx"
#include "vector_io.hxx"

//Takes a vector which should be interpolated
//returns a vector that is the interpolated version of the original vector
//This vector will have a size == input.size() + (input.size()-1)
vector<double> interpolate(const vector<double> &input);


//multiplies the first value with first value, second value with second value, etc. of two arrays (must be same length)
vector<double> mult2arrays(const vector<double> &array1, const vector<double> &array2);

vector<double> Discrete_Convolution_2_Odd(vector<double> list1);

vector<double> Completeness(double Mag_min, double Mag_max, double interval);


vector<double> Completeness_2(vector<double> CC);

//returns chi-squared value
//DOES NOT WORK WHEN THERE ARE ZEROS IN THE EXPECTED (INPUT) DATA
double chi_squared(vector<double> observed, vector<double> expected);

extern vector<double> objective_function_t2;
extern vector<double> storage;
extern double fitness;
double objective_function(const vector<double> &t1);

void optimize(vector<double> t1);

vector<double> startingFitEfficiency(vector<double> startingData);

void recordTranformation(string filePath, vector<double> transformedData, vector<double> startingFit);

void record(string filePath, int line, vector<double> data);


vector<double> getRealTransformedData (int line1, string filePath);
//Takes input data from files with leading zeros
//Remove leading and trailing zeros and return clean data from 16 to 22.5
vector<double> trimData(const vector<double>& input);

void run(int argc, char** argv);

#endif //__DENSITY_TO_STAR_COUNTS_H__
