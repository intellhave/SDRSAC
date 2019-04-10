#include "mex.h"
#include "KDTree.h"

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){   
	// check the arguments
	if( nrhs!=1 || !mxIsNumeric(prhs[0]) )
		mexErrMsgTxt("varargin{1} must be a valid kdtree pointer\n");
	
	// retrieve the tree pointer
    KDTree* tree = KDTree::retrieve_pointer(prhs[0]);
    tree -> ~KDTree();
}
