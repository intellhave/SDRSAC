// this is actually just a single query

#include "KDTree.h"
#include "mex.h"

void retrieve_tree( const mxArray* matptr, KDTree* & tree){
    // retrieve pointer from the MX form
    double* pointer0 = mxGetPr(matptr);
    // check that I actually received something
    if( pointer0 == NULL )
        mexErrMsgTxt("vararg{1} must be a valid k-D tree pointer\n");
    // convert it to "long" datatype (good for addresses)
    intptr_t pointer1 = (intptr_t) pointer0[0];
    // convert it to "KDTree"
    tree = (KDTree*) pointer1;
    // check that I actually received something
    if( tree == NULL )
        mexErrMsgTxt("vararg{1} must be a valid k-D tree pointer\n");
    if( tree -> ndims() <= 0 )
        mexErrMsgTxt("the k-D tree must have k>0"); 
}
void retrieve_point( const mxArray* matptr, vector<double>& point ){
    // check that I actually received something
    if( matptr == NULL )
        mexErrMsgTxt("vararg{2} must be a [kxN] matrix of data\n");
    
    // check dimensions
    if( ( mxGetM(matptr)!=1 || mxGetN(matptr)!=point.size() ) &&
    	( mxGetN(matptr)!=1 || mxGetM(matptr)!=point.size() ) )
    	mexErrMsgTxt("vararg{2} must be a [kx1] or a [1xk] point\n");  
    
    // retrieve point
    double* data = mxGetPr(matptr);
    
    for(int dim=0; dim < point.size(); dim++)
		point[dim] = data[dim];
}
void retrieve_radius( const mxArray* matptr, double& radius ){
    // check that I actually received something
    if( matptr == NULL )
        mexErrMsgTxt("vararg{3} must be a scalar\n");

    if( 1 != mxGetM(matptr) || !mxIsNumeric(matptr) || 1 != mxGetN(matptr) )
    	mexErrMsgTxt("vararg{3} must be a scalar\n");    
    
    // retrieve point
	radius = mxGetScalar(matptr);	
}
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){
	// chec number of arguments
	if( nrhs!=3 )
		mexErrMsgTxt("This function requires 3 arguments\n");
	if( !mxIsNumeric(prhs[0]) )
		mexErrMsgTxt("varargin{0} must be a valid kdtree pointer\n");
	if( !mxIsNumeric(prhs[1]) )
		mexErrMsgTxt("varargin{1} must be a query point\n");
	if( !mxIsNumeric(prhs[2]) )
		mexErrMsgTxt("varargin{2} must be a double (radius)\n");
	
	// retrieve the tree pointer
    KDTree* tree;
    retrieve_tree( prhs[0], tree ); 
    // retrieve the query point
    vector<double> point(tree->ndims(),0);
    retrieve_point(prhs[1], point);
    // retrieve the radius
    double radius = 0;
    retrieve_radius(prhs[2], radius);
    
    vector<int> idxsInRange;
    vector<double> dists;
    tree->ball_query( point, radius, idxsInRange, dists );
    
    if( nlhs == 1 ){
	    // convert vector<int> back in matlab format
	    plhs[0] = mxCreateDoubleMatrix(idxsInRange.size(), 1, mxREAL);
	    double* indexes = mxGetPr(plhs[0]);
	    for (int i=0; i < idxsInRange.size(); i++)
	    	indexes[ i ] = idxsInRange[i] + 1;
    } else 
    if( nlhs == 2 ){
	    // convert vector<int> back in matlab format
	    plhs[0] = mxCreateDoubleMatrix(idxsInRange.size(), 1, mxREAL);
	    double* indexes = mxGetPr(plhs[0]);
	    for (int i=0; i < idxsInRange.size(); i++)
	    	indexes[ i ] = idxsInRange[i] + 1;
	    
	    // convert vector<double> back in matlab format
	    plhs[1] = mxCreateDoubleMatrix(dists.size(), 1, mxREAL);
	    double* distances = mxGetPr(plhs[1]);
	    for (int i=0; i < dists.size(); i++)
	    	distances[ i ] = dists[i];
    	
    }
    else
    	mexErrMsgTxt("provide either one or two output parameters.");
}
