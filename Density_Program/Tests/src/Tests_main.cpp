#include "Density_toStarCounts_Tests.h"
using namespace std;

//Integration Tests
//
void testInputvsResults(const vector<double> &t1)
{
	//vector<double> starcounts = interpolate(t1);
		
	vector<double> starcounts = t1;
	//Calculate magnitude ranges
	double lowerbound = 16;
	double upperbound = 16 + double(starcounts.size())/2.0;
	
    vector<double> convolved = Discrete_Convolution_2_Odd(starcounts);
    
	vector<double> CC = Completeness(lowerbound, upperbound, .5);

	//because the convolution adds extra bins on the end, this removes the bins
	//to make the two vectors have the same size
	while (convolved.size() != CC.size())
	{
		convolved.pop_back();
		convolved.erase(convolved.begin());
	}
	
	vector<double> final_star_count = mult2arrays(convolved, CC);
    
    cout << "Input: ";
    for(size_t i = 0; i < t1.size(); i++)
    {
        cout << t1[i] << ", ";
    }
    cout << endl;
    cout << "Processed: ";
    for(size_t i = 0; i < final_star_count.size(); i++)
    {
        cout << final_star_count[i] << ", ";
    }
    cout << endl;
    
}

void testDelimiterFinding()
{
    //Test for "unsigned int findLastDelimiter(string input, char delimiter)"
    cout << "\nTest findLastDelimiter\n====================================\n\n";
    unsigned int failed = 0;
    failed += (findLastDelimiter("../../Work/12/folderName", '/') != 13);
    failed += (findLastDelimiter("../../Work/12/folderName/", '/') != 24);
    failed += (findLastDelimiter("/", '/') != 0);
    failed += (findLastDelimiter("", '/') != 0);
    failed += (findLastDelimiter("..*..*Work*12*folderName", '/') != 24);
    failed += (findLastDelimiter("/.*..*Work*12*folderName", '/') != 0);
    if ( failed )
    {
        cout << "findLastDelimiter() failed " << failed << "/6 tests." << endl;
    }
    else
    {
        cout << "findLastDelimiter() passed all tests." << endl;
    }
}

void testInterpolate()
{
	//Testing interpolate
	vector<double> t1 {0,1,2,3,4,5,6,7,8,9,10};
	vector<double> r1 = interpolate(t1); 
	std::cout<< "\nTest Interpolate\n====================================\n\n";
	
	std::cout<< "Original Before Interpolation; Size: " <<t1.size()<< std::endl;
	for (size_t i = 0; i < t1.size(); i++)
	{
		std::cout<<t1[i]<<", ";
	}
	std::cout <<"\nAfter Interpolation; Size: "<< r1.size() <<std::endl;	
	for (size_t i = 0; i < r1.size(); i++) 
	{
		std::cout<< r1[i] << ", "; 
	}
	std::cout << "\n";
	
	int fail = 0; 
	if (r1.size() != 2 * t1.size() - 1)
		std::cout<<"The Interpolation function does not apply the correct change in size."<<std::endl;
	for (size_t i = 0; i < r1.size(); i++)
	{
		if (i%2 == 0)
		{
			if (r1[i] != t1[i/2])
				fail = 1;
		}
		else 
		{
			double avg = t1[i/2] + t1[i/2 + 1];
			avg /= 2;
			if (r1[i] != avg)
				fail = 1;
		}
	}
	if (fail == 1)
		std::cout<<"The interpolation function does not apply the correct averaging in between vectors."<<std::endl;
	else
		std::cout<<"testInterpolate passed all tests"<<std::endl;
	
}

void printArray(const vector<double>& ex)
{
	for (size_t i = 0; i < ex.size(); i++)
	{
		cout<<ex[i];
		if (i != ex.size()-1)
			cout<<", ";
		else
			cout<<"\n";
	}
}

void testMult2arrays()
{
	std::cout<<"\nTest mult2arrays\n====================================\n\n";
	vector<double> a1 = {123.1234, 9, 0, .000193, -12314.49, 84372.1231};
	vector<double> a2 = {1, 92021, .5, .1131, 13, -19392};
	vector<double> diff = {0};
	
	//if the size is different, return the first array size, with all elements 
	//equal to 0
	vector<double> res = mult2arrays(a2, diff);
	bool flag1 = false; 
	//checks the returned result has the same size as the first array
	flag1 = !(res.size() == a2.size()); 
	for(size_t i = 0; i < res.size(); i++)
	{
		//checks that the result is the same as the first array
		if (res[i] != 0)
			flag1 = true;
	}
	if (flag1)
		std::cout<<"Wrong output when the two arrays are different"<<std::endl;
	
	bool flag2 = false;
	
	cout<<"Array 1: ";
	printArray(a1);
	cout<<"Array 2: ";
	printArray(a2);

	//if the size is the same, return a_i * b_i for all i
	res = mult2arrays(a1, a2);
	
	cout<<"Result Array: ";
	printArray(res);
	for (size_t i = 0; i < res.size(); i++)
	{
		//checks each element is the multiplied of the other two elements
		if (res[i] != a1[i] * a2[i])
			flag2 = true; 
	}
	if (flag2)
		std::cout<<"Mult2arrays has the wrong behavior"<<std::endl;
	if (flag1 or flag2)
		std::cout<<"testMult2arrays has FAILED"<<std::endl; 
	else
		std::cout<<"testMult2arrays passed all tests"<<std::endl;
}

bool testDiscrete_Convolution_2_Helper(const vector<double>& input)
{
	double presum = 0;
	std::cout<<"Input Vector: ";
	for(size_t i = 0; i < input.size(); i++)
	{
		presum += input[i];
		std::cout << input[i];
		if (i < input.size() - 1)
			std::cout<<",";
	}
	vector<double> temp = Discrete_Convolution_2_Odd(input); 
	std::cout<<"\nOutput Vector: ";
	double postsum = 0;
	for(size_t i = 0; i < temp.size(); i++)
	{
		std::cout << temp[i];
		postsum += temp[i];
		if (i < temp.size() - 1)
			std::cout<<",";
	}
	std::cout << "\nPreSum: "<<presum<<", PostSum: "<< postsum << "\n\n";
	if (input.size() != temp.size())
		return true;
	else
		return false;
}	

void testDiscrete_Convolution_2_Odd()
{
	std::cout << "\nTest Discrete_Convolution_2_Odd\n====================================\n\n";
	vector<double> temp;
	vector<double> input = {0,0,0,0,0,0,0,0,0,0,0,0,0};
	bool flagsize = false; 
	for(size_t val = 0; val < input.size(); val++)
	{
		input[val] = 1;
		if (val > 0)
			input[val-1]=0;
		flagsize = testDiscrete_Convolution_2_Helper(input) || flagsize;
	}
	
	input = {0,3,4,1,0,2,6,4,4,0,0,3,0};
	flagsize = testDiscrete_Convolution_2_Helper(input);
	std::cout<<"Checking whether gaussians are normalized.\n\nGaussian Areas:"<<std::endl; 
	//tests whether gaussian function is normalized
	vector< vector <double> > vectorGaussian = {{0.000146222,0.00152518,0.0104476,0.0470463,0.139381,0.271853,0.343851,0.169326,0.0163913,0.000275136,0.000000749684},{0.000146292,0.00152567,0.0104498,0.0470513,0.139386,0.271852,0.343844,0.169322,0.016391,0.00027513,0.000000749667},{0.00014647,0.00152694,0.0104552,0.0470641,0.139399,0.271849,0.343826,0.169312,0.01639,0.000275115,0.000000749625},{0.000146993,0.00153067,0.0104711,0.0471017,0.139438,0.271838,0.34377,0.169284,0.0163873,0.000275069,0.000000749499},{0.000148947,0.00154456,0.0105301,0.0472411,0.139583,0.271799,0.343566,0.169179,0.0163771,0.000274897,0.000000749033},{0.000158861,0.00161417,0.0108227,0.0479266,0.140287,0.271603,0.342566,0.168663,0.0163271,0.000274059,0.000000746749},{0.000240385,0.0021425,0.0129031,0.0525479,0.144804,0.270141,0.335936,0.165252,0.015997,0.000268518,0.000000731649},{0.00202997,0.00913824,0.0314354,0.0826496,0.166109,0.255227,0.294742,0.144323,0.013971,0.00023451,0.000000638985},{0.0124662,0.0306631,0.0642212,0.114533,0.173935,0.224931,0.243223,0.118616,0.0114824,0.000192738,0.000000525168},{0.0145973,0.0339889,0.0680515,0.117162,0.173455,0.220825,0.237361,0.115715,0.0112016,0.000188025,0.000000512324},{0.0146547,0.0340756,0.0681486,0.117226,0.17344,0.220719,0.237212,0.115641,0.0111945,0.000187905,0.000000511997},{0.0146552,0.0340765,0.0681495,0.117226,0.173439,0.220718,0.23721,0.11564,0.0111944,0.000187903,0.000000511994},{0.0146552,0.0340765,0.0681495,0.117226,0.173439,0.220718,0.23721,0.11564,0.0111944,0.000187903,0.000000511994}};
    for(size_t i = 0; i < vectorGaussian.size(); i++)
	{
		double total = 0;
		for(size_t j = 0; j < vectorGaussian[i].size(); j++)
		{
			total += vectorGaussian[i][j];
		}
		if (i != vectorGaussian.size()-1)
			std::cout<<total<< ", ";
	}
    std::cout << std::endl;
	if (flagsize)
		std::cout<<"The Discrete_Convolution_2_Odd() returns a bins of the that is of different size to that of the input?!?!"<<std::endl;
	cout<<"Finished all tests for Discrete_Convolution_2_Odd\n";
}

void testCompleteness()
{
	std::cout << "\nTest Completeness\n====================================\n\n";

    vector <double> temp = Completeness(16, (22.5), .5);
    cout<<"Completeness from mag 16 to 22 in .5 steps"<<endl;
    printArray(temp);
    std::cout << std::endl;
}


void unitTests()
{
	testInterpolate();
    testMult2arrays();
	testDiscrete_Convolution_2_Odd();
	testCompleteness();
	


    testDelimiterFinding();

}

//I don't think I actually want to use this.
bool compareVectors(vector<double> vec1, vector<double> vec2){
    if(vec1.size() != vec2.size())
    {
        return false; 
    }
    
    //Same size is guaranteed ^^^^
    for(size_t i = 0; i < vec1.size(); i++)
    {
    
        if(vec1[i] != vec2[i])
        {
            return false;
        }
    }
    return true;
}

//Unit Tests

void runTests(int argc, char** argv) 
{
    unitTests();
    
    testInputvsResults({1105.97,1244.78,1383.58,2508.33,3633.07,4097.92,4562.76,21408.,38253.3,60745.9,83238.6,51836.7,20434.7});
    
    return;
}
