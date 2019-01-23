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

#include "tao/asynchronous_algorithms/particle_swarm.hxx"
#include "tao/asynchronous_algorithms/differential_evolution.hxx"

#include "tao/synchronous_algorithms/parameter_sweep.hxx"
#include "tao/synchronous_algorithms/synchronous_gradient_descent.hxx"
#include "tao/synchronous_algorithms/synchronous_newton_method.hxx"

//from undvc_common
#include "undvc_common/arguments.hxx"
#include "undvc_common/vector_io.hxx"
using namespace std;


//Takes a vector which should be interpolated
//returns a vector that is the interpolated version of the original vector
//This vector will have a size == input.size() + (input.size()-1)
vector<double> interpolate(const vector<double> &input) 
{
    vector<double> result = vector<double>(input.size()+(input.size()-1), 0.0);
	//Insert input values into result in every other position so we can interpolate
    for (int i = 0; i < result.size(); i += 2)
    {
        result[i] = input[i/2];
    }
    //Interpolate the values from input using linear interpolation
    for (int i = 1; i < result.size(); i += 2)
    {
        result[i] = (result[i-1]+result[i+1])/2.0;
    }
    
    return result;
}


//multiplies the first value with first value, second value with second value, etc. of two arrays (must be same length)
vector<double> mult2arrays(const vector<double> &array1, const vector<double> &array2) //ampersands are the address of a variable, which is why they are constants here (i think)
{
	int size1 = array1.size();
	int size2 = array2.size();
    
	vector<double> result (size1, 0);
    
	if (size1 != size2)
	{
		return result;
	}
	
	for (int i = 0; i < size1; i++)
	{
		result[i] = array1[i] * array2[i];
	}
	
	return result; 
}

vector<double> Discrete_Convolution_2_Odd(vector<double> list1)
{
    vector<double> histogram1 = list1;

    vector< vector <double> > vectorGaussian = {{0.000146222,0.00152518,0.0104476,0.0470463,0.139381,0.271853,0.343851,0.169326,0.0163913,0.000275136,0.000000749684},{0.000146292,0.00152567,0.0104498,0.0470513,0.139386,0.271852,0.343844,0.169322,0.016391,0.00027513,0.000000749667},{0.00014647,0.00152694,0.0104552,0.0470641,0.139399,0.271849,0.343826,0.169312,0.01639,0.000275115,0.000000749625},{0.000146993,0.00153067,0.0104711,0.0471017,0.139438,0.271838,0.34377,0.169284,0.0163873,0.000275069,0.000000749499},{0.000148947,0.00154456,0.0105301,0.0472411,0.139583,0.271799,0.343566,0.169179,0.0163771,0.000274897,0.000000749033},{0.000158861,0.00161417,0.0108227,0.0479266,0.140287,0.271603,0.342566,0.168663,0.0163271,0.000274059,0.000000746749},{0.000240385,0.0021425,0.0129031,0.0525479,0.144804,0.270141,0.335936,0.165252,0.015997,0.000268518,0.000000731649},{0.00202997,0.00913824,0.0314354,0.0826496,0.166109,0.255227,0.294742,0.144323,0.013971,0.00023451,0.000000638985},{0.0124662,0.0306631,0.0642212,0.114533,0.173935,0.224931,0.243223,0.118616,0.0114824,0.000192738,0.000000525168},{0.0145973,0.0339889,0.0680515,0.117162,0.173455,0.220825,0.237361,0.115715,0.0112016,0.000188025,0.000000512324},{0.0146547,0.0340756,0.0681486,0.117226,0.17344,0.220719,0.237212,0.115641,0.0111945,0.000187905,0.000000511997},{0.0146552,0.0340765,0.0681495,0.117226,0.173439,0.220718,0.23721,0.11564,0.0111944,0.000187903,0.000000511994},{0.0146552,0.0340765,0.0681495,0.117226,0.173439,0.220718,0.23721,0.11564,0.0111944,0.000187903,0.000000511994}};
    
    vector<double> result;
    
    vector<bool> included(vectorGaussian[0].size(), true);
    
    int edges = (int)(vectorGaussian[0].size()/2);
    
    for (unsigned int i = 0; i < histogram1.size(); i++) //changed from unsigned int to int
    {
        for (unsigned int t = 0; t < vectorGaussian[i].size(); t++)
            included[t] = true;
        
        //checks if the gaussian bars lie over the stars
        for (int j = -(edges); j <= edges; j++)
        {
            if ((i+j < 0) || (i+j >= (int) histogram1.size()))
            {
                //getting a warning because an expression with an unsigned variable in it will always be greater than zero; so this statment only depends on the second option
                // when run by itself, when i = 0 and j = -2, i + j = 4294967294, which is not -2
                included[j+edges] = false; //changed from 2
            }
        }
        
        double overall = 1.;
        
        for (int m = 0; m < vectorGaussian[i].size(); m++)
        {
            if (!included[m])
            {
                overall = 0.;
            }
        }
        
        //finds the amount of area under the remaining bars
        if (overall < .001)
        {
            overall = 0.;
            for (unsigned int alpha = 0; alpha < vectorGaussian[i].size(); alpha++)
            {
                overall += ((double) included[alpha]) * vectorGaussian[i][alpha];
            }
        }
        double total = 0.;
        
        
        //finds total to add to convolution, when gaussian is not over the stars does not count anything.
        //check this for edges
        for(unsigned int k = 0; k < vectorGaussian[i].size(); k++)
        {
            if (!included[k])
                continue;
            total += histogram1[i+k-edges] * vectorGaussian[i][k] / overall;
        }
        
        //adds total to result list of values
        result.push_back(total);
    }
    
    return result;
}


vector<double> Completeness(double Mag_min, double Mag_max, double interval)
{
	double s0 = .9402;
	double s1 = 1.6171;
	double s2 = 23.5877;
	int iterations = (int)((Mag_max - Mag_min) / interval);
	
	vector<double> CC;
    
    
    //adds values of completeness for each iteration in the array of stars
	for (double i = 0; i < iterations; i++)
	{
		double value = s0 / (exp(s1*(i * interval + Mag_min + (interval/2.0) - s2)) + 1); // gives 16.25 - s2, 16.75 - s2, 17.25 - s2, etc. it is using the midpointof each bin as magnitude
		CC.push_back(value);
	}
	return CC; 
}


vector<double> Completeness_2(vector<double> CC)
{
    //Initialize the array to the completeness coefficient for each bin
    vector<double> vectorCompleteness = {0.986088,0.972016,0.957948,0.944999,0.934449,0.927228,0.922852,0.917635,0.902512,0.861945,0.777409,0.640076,0.471073};
    assert(CC.size() == vectorCompleteness.size());
    
    //Element-wise multiply the input array by the completeness
    for (int i = 0; i < vectorCompleteness.size(); i++)
    {
        vectorCompleteness[i] *= CC[i];  
    }

    return vectorCompleteness;
    
}

//returns chi-squared value
//DOES NOT WORK WHEN THERE ARE ZEROS IN THE EXPECTED (INPUT) DATA
double chi_squared(vector<double> observed, vector<double> expected)
{
	if (observed.size() != expected.size())
		return -1;
    
	double chi_sq = 0;
    
	for (unsigned int i = 0; i < observed.size(); i++)
	{
		chi_sq += ((observed[i] - expected[i])*(observed[i] - expected[i]))/expected[i];
	}
	return chi_sq; 
}


//function that takes input (2 histograms with lengths) with output (chi-squared value)
//t1 is the guess of star counts, t2 is the expected star counts
//Temporarily add global variable since TAO can only handle objective functions with 1 argument
//t1 should be of a length that when interpolated it is the same length as t2
vector<double> objective_function_t2;
vector<double> storage;
double fitness;

double objective_function(const vector<double> &t1)
{	
	//Interpolate guess and store for output
	vector<double> starcounts = interpolate(t1);
    storage = starcounts;
		
	//Calculate magnitude ranges
	double lowerbound = 16;
	double upperbound = 16 + double(starcounts.size())/2.0;
	
    vector<double> convolved = Discrete_Convolution_2_Odd(starcounts);
    
	vector<double> CC = Completeness_2(Completeness(lowerbound, upperbound, .5));

	//because the convolution adds extra bins on the end, this removes the bins
	//to make the two vectors have the same size
	while (convolved.size() != CC.size())
	{
		convolved.pop_back();
		convolved.erase(convolved.begin());
	}
	
	vector<double> final_star_count = mult2arrays(convolved, CC);

	double result = chi_squared(objective_function_t2, final_star_count);
    fitness = result;

	return -(result);
}

void optimize(vector<double> t1)
{
    
    std::vector<double> min_bound(t1.size(), 0.5); //recall this makes a vector of size t1.size() full of 0.0
    std::vector<double> max_bound(t1.size(), 100000.0); // was a million
    
    
    vector<double> step_size(t1.size(), 1); //would decreasing the step size help the smoothness problem? (step size was originally 1)
    vector<double> starting_point = t1;
    
	vector<double> final_parameters;
    
	double final_fitness = 0;
	vector <string> args;
	args.push_back("--min_improvement");
	
	
	args.push_back("1e-10");//from e-10
	args.push_back("--gd_quiet");//from e-10
	args.push_back("--max_iterations");
	
	
	args.push_back("1000000");
    synchronous_gradient_descent(args, objective_function, min_bound, max_bound, starting_point, step_size, final_parameters, final_fitness); //added min_bound and max_bound
}

vector<double> startingFitEfficiency(vector<double> startingData)
{
    vector<double> startingFit;
    vector<double> CC = Completeness_2(Completeness(16., (startingData.size()/2.0 + 16.), .5));
    
    for (int i = 0; i < startingData.size(); i++)
    {
        startingFit.push_back(startingData[i]/CC[i]);
    }
    
    return startingFit;
}

void recordTranformation(string filePath, vector<double> transformedData, vector<double> startingFit)
{
    ofstream recorded;
    recorded.open(filePath, ofstream::out | ofstream::trunc);
    
    recorded << "Number," << "16," << "16.5," << "17," << "17.5," << "18," << "18.5," << "19," << "19.5," << "20," << "20.5," << "21," << "21.5," << "22," << "22.5" << endl;
    
    recorded << "Transformed Data:";
    for (unsigned int j = 0; j < transformedData.size(); j++)
    {
        recorded << "," << transformedData[j];
    }
    
    recorded << endl << "Starting Fit:";
    for (unsigned int j = 0; j < startingFit.size(); j++)
    {
        recorded << "," << startingFit[j];
    }
    
    recorded.close();
}

void record(string filePath, int line, vector<double> data) //had vector<double> data
{
    ofstream recorded;
    recorded.open(filePath, ofstream::out | ofstream::app);
    
    recorded << endl << to_string(line);
    for (unsigned int j = 0; j < data.size(); j++)
    {
        recorded << "," << data[j];
    }
    
    recorded << "," << fitness;
    
    recorded.close();
}

//Read data from file
vector<double> getRealTransformedData (int line1, string filePath)
{
    
    ifstream infile(filePath); // for example
    string line = "";
    vector<double> Attempt1;
    
    if (infile.is_open())
    {
        while (getline(infile, line, '\n'))
        {
            stringstream strstr(line);
            string word = "";
            
            getline(strstr,word, ',');
            if (word == to_string(line1))
            {
                while (getline(strstr,word, ','))
                {
                    Attempt1.push_back(stod(word));
                    
                }
            }
            
        }
    }
    else
    {
        cout << "Failed to open file" << endl;
    }
    

/*    
    cout << "number (" << Attempt1.size() << "): ";
    for (int i = 0; i < Attempt1.size(); i++)
    {
    cout << Attempt1[i] << ", ";
    }
    
    cout << endl;
*/
    
    return Attempt1;

}


void testInputvsResults()
{
	vector<double> starcounts = interpolate(t1);
    storage = starcounts;
		
	//Calculate magnitude ranges
	double lowerbound = 16;
	double upperbound = 16 + double(starcounts.size())/2.0;
	
    vector<double> convolved = Discrete_Convolution_2_Odd(starcounts);
    
	vector<double> CC = Completeness_2(Completeness(lowerbound, upperbound, .5));

	//because the convolution adds extra bins on the end, this removes the bins
	//to make the two vectors have the same size
	while (convolved.size() != CC.size())
	{
		convolved.pop_back();
		convolved.erase(convolved.begin());
	}
	
	vector<double> final_star_count = mult2arrays(convolved, CC);

}



void runTests() 
{

    vector<double> temp = Completeness_2(vector<double>(13, 1.0));
    std::cout << "Test Completeness_2" << std::endl;
    for(size_t i = 0; i < temp.size(); i++)
    {
        std::cout << temp[i] << ", ";
    }
    std::cout << std::endl;
    
    temp = Completeness(16, (22.5), .5);
    std::cout << "Test Completeness" << std::endl;
    for(size_t i = 0; i < temp.size(); i++)
    {
        std::cout << temp[i] << ", ";
    }
    std::cout << std::endl;
    
    temp = Discrete_Convolution_2_Odd({0,0,0,0,0,0,1,1,0,0,0,0,0});
    std::cout << "Test Discrete_Convolution_2_Odd" << std::endl;
    for(size_t i = 0; i < temp.size(); i++)
    {
        std::cout << temp[i] << ", ";
    }
    std::cout << std::endl;
    
    
    
    
    
    return;
}


//Takes input data from files with leading zeros
//Remove leading and trailing zeros and return clean data from 16 to 22.5
vector<double> trimData(const vector<double>& input)
{
    vector<double> result = vector<double>(13, 0.0);
    for (int i = 0; i < result.size(); i++)
    {
        cout << input[i+32] << ", ";
        result[i] = input[i+32];
    }
    cout << endl;

    return result;
}


int main(int argc, char* argv[])
{
    //runTests();
    for (int i = 2; i < argc; i++)
    {
        string filePathInput = argv[i];
        string filePathOutput = string("./JakeWork/PannStars_Results/") + &argv[i][26];
        //string filePathStartingFit = argv[1];
        
        //Input data must have odd number of bins
        vector<double> realTransformedData = getRealTransformedData(0, filePathInput);
        if((argv[1]) == string("1"))
        {
            realTransformedData = trimData(realTransformedData);
        }
        assert(realTransformedData.size()%2);
        
        unsigned int numberOfFitParameters = realTransformedData.size()-(realTransformedData.size()/2); 
        vector<double> startingFit = vector<double>(numberOfFitParameters, 0.0);
        //vector<double> efficiencyCorrected = startingFitEfficiency(realTransformedData);

        int totalStars = 0;
        for(unsigned int j = 0; j < startingFit.size(); j++)
        {
            //startingFit[j] = efficiencyCorrected[j*2];
            startingFit[j] = realTransformedData[j*2];
            totalStars += startingFit[j];
        }
        
        if(totalStars > 0)
        {
            std::cerr << "Input: " << filePathInput << std::endl;
            std::cerr << "Output: " << filePathOutput << std::endl;
            std::cerr << "Input Number: " << i << std::endl;
            //recordTranformation(filePathStartingFit, realTransformedData, startingFit); //got rid of trial
            
            objective_function_t2 = realTransformedData;
            
            optimize(startingFit);
            
            //Clear old file
            ofstream recorded;
            recorded.open(filePathOutput, ofstream::out);
            recorded.close();
            
            record(filePathOutput, 0, storage); //got rid of trial
            record(filePathOutput, 1, startingFit); //got rid of trial
            record(filePathOutput, 2, realTransformedData); //got rid of trial
        }
        else
        {
            //cout << "No stars in this histogram" << endl;
        }    
    }
    return 0;
}
