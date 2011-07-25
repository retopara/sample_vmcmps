#include <iostream>

using namespace std;

int flip(int spin)
{
	int flip;
	if(spin == 0) flip = 1;
	if(spin == 1) flip = 0;
	return flip;
}
