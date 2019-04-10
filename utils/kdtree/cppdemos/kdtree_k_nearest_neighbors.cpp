#include "KDTree.h"

int test1(){
	int N = 30;
	vector< Point > A(N, vector<double>(2,0));
	A[0][0] = 0.51; A[0][1] = 0.79;
	A[1][0] = 0.46; A[1][1] = 0.41;
	A[2][0] = 0.35; A[2][1] = 0.12;
	A[3][0] = 0.10; A[3][1] = 0.63;
	A[4][0] = 0.43; A[4][1] = 0.86;
	A[5][0] = 0.71; A[5][1] = 0.16;
	A[6][0] = 0.12; A[6][1] = 0.60;
	A[7][0] = 0.08; A[7][1] = 0.12;
	A[8][0] = 0.37; A[8][1] = 0.63;
	A[9][0] = 0.03; A[9][1] = 0.84;
	A[10][0] = 0.19; A[10][1] = 0.94;
	A[11][0] = 0.47; A[11][1] = 0.42;
	A[12][0] = 0.14; A[12][1] = 0.27;
	A[13][0] = 0.72; A[13][1] = 0.93;
	A[14][0] = 0.66; A[14][1] = 0.92;
	A[15][0] = 0.43; A[15][1] = 0.54;
	A[16][0] = 0.45; A[16][1] = 0.81;
	A[17][0] = 0.51; A[17][1] = 0.17;
	A[18][0] = 0.53; A[18][1] = 0.32;
	A[19][0] = 0.57; A[19][1] = 0.66;
	A[20][0] = 0.36; A[20][1] = 0.00;
	A[21][0] = 0.34; A[21][1] = 0.63;
	A[22][0] = 0.17; A[22][1] = 0.79;
	A[23][0] = 0.09; A[23][1] = 0.29;
	A[24][0] = 0.39; A[24][1] = 0.79;
	A[25][0] = 0.80; A[25][1] = 0.22;
	A[26][0] = 0.01; A[26][1] = 0.40;
	A[27][0] = 0.23; A[27][1] = 0.80;
	A[28][0] = 0.93; A[28][1] = 0.86;
	A[29][0] = 0.23; A[29][1] = 0.14;

	KDTree* tree = new KDTree( A );
	cout << "tree created" << endl;
	tree -> print_tree();
	// create query
	vector<double> point(2,0);
//	point[0]=.45; point[1]=.8;
	point[0]=.5; point[1]=.5;


	// kNN query
	vector<int> idxsInRange;
    vector<double> dists;
    tree->k_closest_points( point, 3, idxsInRange, dists );
    cout << "query data executed" << endl;

    cout << "points in range: ";
    for (unsigned int i=0; i < idxsInRange.size(); i++)
		cout << idxsInRange[i] + 1 << " ";
    cout << endl;

	return 0;
}
int test2(){
	int N = 10;
	vector< Point > A(N, vector<double>(2,0));
	A[0][0] = 0.10; A[0][1] = 0.54;
	A[1][0] = 0.68; A[1][1] = 0.77;
	A[2][0] = 0.80; A[2][1] = 0.16;
	A[3][0] = 0.76; A[3][1] = 0.74;
	A[4][0] = 0.47; A[4][1] = 0.29;
	A[5][0] = 0.67; A[5][1] = 0.45;
	A[6][0] = 0.93; A[6][1] = 0.57;
	A[7][0] = 0.62; A[7][1] = 0.07;
	A[8][0] = 0.95; A[8][1] = 0.22;
	A[9][0] = 0.27; A[9][1] = 0.58;

	KDTree* tree = new KDTree( A );
	tree -> print_tree();
	vector<double> point(2,0);
	point[0]=.8; point[1]=.1;
	int k = 1;

	// kNN query
	vector<int> idxsInRange;
    vector<double> dists;
    tree->k_closest_points( point, k, idxsInRange, dists );

    cout << "points in range: ";
    for (unsigned int i=0; i < idxsInRange.size(); i++)
		cout << idxsInRange[i] + 1 << " ";
    cout << endl;

	return 0;
}

int main(){
	// test1();
	test2();
}
