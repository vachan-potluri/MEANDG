/** @file */
#include "sysInclude.h"
#include "Tests.h"
using namespace std;


int main(){

	srand(int(time(0)));

	//runAllTests();
	Test::TestPoint();
	Test::TestTensor();

	return 0;
}
