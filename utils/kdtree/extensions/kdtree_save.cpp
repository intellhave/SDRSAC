#include "KDTree.h"

#ifndef CPPONLY
#include "mex.h"

// This function code was taken from "kdtree_nearest_neighbor.cpp":
void retrieve_tree( const mxArray* matptr, KDTree* & tree)
{
    // retrieve pointer from the MX form
    double* pointer0 = mxGetPr(matptr);
    // check that I actually received something
    if( pointer0 == NULL )
        mexErrMsgTxt("vararg{1} must be a valid k-D tree pointer\n");
    // convert it to "long" datatype (good for addresses)
    long pointer1 = (long) pointer0[0];
    // convert it to "KDTree"
    tree = (KDTree*) pointer1;
    // check that I actually received something
    if( tree == NULL )
        mexErrMsgTxt("vararg{1} must be a valid k-D tree pointer\n");
    if( tree -> ndims() <= 0 )
        mexErrMsgTxt("the k-D tree must have k>0"); 
}

// This function was created by Pablo Sala on November 25, 2008
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    if(nlhs != 0 || nrhs != 2 ||
       !mxIsNumeric(prhs[0])  ||  // tree     is an integer
       !mxIsChar(prhs[1]))        // filename is a  string
    {
        mexErrMsgTxt("Incorrect Sintax\n\t kdtree_save(tree, filename),\n\twhere 'tree' is a kd-tree and 'filename' is a string\n");
    }

    // Retrieve the pointer to the kd-tree:
    KDTree* tree;
    retrieve_tree(prhs[0], tree);

    // Retrieve the filename:
    const char* filename = mxArrayToString(prhs[1]);

    // Save the kd-tree:
    tree->save(filename);
}

#endif
