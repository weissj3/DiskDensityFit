#include "Density_to_Star_Counts.h"
#ifdef TEST
#include "Density_toStarCounts_Tests.h"
#endif //TEST

using namespace std;

int main(int argc, char* argv[])
{
    #ifdef TEST
    runTests(argc, argv);
    #else
    run(argc, argv);
    #endif //TEST

    
    return 0;
}
