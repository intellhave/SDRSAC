#include "KDTree.h"
#include "mex.h"
#include "matrix.h" //isNaN/isinf

// matlab entry point
void retrieve_data( const mxArray* matptr, vector< vector<double> >& dataV, int& npoints, int& ndims){
    // retrieve pointer from the MX form
    double* data = mxGetPr(matptr);
    // check that I actually received something
    if( data == NULL )
        mexErrMsgTxt("vararg{2} must be a [kxN] matrix of data\n");
    
    // retrieve amount of points
    npoints = mxGetM(matptr);
    ndims   = mxGetN(matptr);

    // Make sure !nan & !inf
    for( int i=0; i<npoints*ndims; i++ ){
        if( mxIsNaN( data[i] ) ) mexErrMsgTxt("input data contains NAN values.");
        if( mxIsInf( data[i] ) ) mexErrMsgTxt("input data contains INF values!");
    }
    
    // FILL THE DATA STRUCTURES
	dataV.resize(npoints, vector<double>(ndims));
	for( int i=0; i<npoints; i++ )
        for( int j=0; j<ndims; j++ )
			dataV[i][j] = data[ i + j*npoints ]; 
}
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){   
	// check input
	if( nrhs != 1 || !mxIsNumeric(prhs[0]) )
		mexErrMsgTxt("A unique [kxN] matrix of points should be passed.\n");
	   
    // retrieve the data
    vector< vector<double> > input_data;
    int npoints;
    int ndims;
    retrieve_data( prhs[0], input_data, npoints, ndims );
    // printf("npoints %d ndims %d\n", npoints, ndims);
    
    // fill the k-D tree
	KDTree* tree = new KDTree( input_data );	
	
	// DEBUG
 	//mexPrintf("npoint %d dimensions %d\n", (int)input_data.size(), (int)input_data[0].size());
	
	// DEBUG
	//for (unsigned int i=0; i < input_data.size(); i++){
	//	for (unsigned int j=0; j < input_data[i].size(); j++)
	//		mexPrintf("%.2f ", input_data[i][j] );
	//	mexPrintf("\n");
	//}

    // return the program a pointer to the created tree
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    double* pointer_to_tree = mxGetPr(plhs[0]);
    pointer_to_tree[0] = (intptr_t) tree;
}
