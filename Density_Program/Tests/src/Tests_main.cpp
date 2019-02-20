#include "Density_toStarCounts_Tests.h"
using namespace std;

//Integration Tests

void testInputvsResults(const vector<double> &t1)
{
	//vector<double> starcounts = interpolate(t1);
		
	vector<double> starcounts = t1;
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
    
    cout << "Input: ";
    for(int i = 0; i < t1.size(); i++)
    {
        cout << t1[i] << ", ";
    }
    cout << endl;
    cout << "Processed: ";
    for(int i = 0; i < final_star_count.size(); i++)
    {
        cout << final_star_count[i] << ", ";
    }
    cout << endl;
    
}

void testDelimiterFinding()
{
    //Test for "unsigned int findLastDelimiter(string input, char delimiter)"
    cout << "Test findLastDelimiter" << endl;
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

void unitTests()
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
    
    testDelimiterFinding();

}

//I don't think I actually want to use this.
bool compareVectors(vector<double> vec1, vector<double> vec2){
    if(vec1.size() != vec2.size())
    {
        return false; 
    }
    
    //Same size is guaranteed ^^^^
    for(int i = 0; i < vec1.size(); i++)
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
