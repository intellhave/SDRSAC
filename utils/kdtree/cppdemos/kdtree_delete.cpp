#include <iostream>
#include "KDTree.h"

using namespace std;

int test1(){
	int N = 8;
	vector< Point > A(N, vector<double>(3,0));
	A[0][0] = 200; A[0][1] = 300; A[0][2] = 400;
	A[1][0] = 100; A[1][1] = 200; A[1][2] = 300;
	A[2][0] = 600; A[2][1] = 100; A[2][2] = 500;
	A[3][0] = 100; A[3][1] = 100; A[3][2] = 300;
	A[4][0] = 700; A[4][1] = 50;  A[4][2] = 300;
	A[5][0] = 700; A[5][1] = 110; A[5][2] = 500;
	A[6][0] = 699; A[6][1] = 51;  A[6][2] = 301;
	A[7][0] = 701; A[7][1] = 50;  A[7][2] = 305;
	KDTree* tree = new KDTree( A );
	tree -> print_tree();
	tree -> ~KDTree();
	cout << "terminated correctly" << endl;
	return 0;
}
int main (int argc, char * const argv[]) {
	return test1(); // simple destruction
}
