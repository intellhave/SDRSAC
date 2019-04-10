#include "KDTree.h"

#ifndef CPPONLY
#include "mex.h"

// Created by Pablo Sala, November 25, 2008

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    if(nlhs != 1 || nrhs != 1 ||
       !mxIsChar(prhs[0]))  // filename is a string
    {
        mexErrMsgTxt("Incorrect Sintax\n\t [tree] = kdtree_load(filename),\n\twhere 'filename' is a string\n");
    }

    // Retrieve the filename:
    const char* filename = mxArrayToString(prhs[0]);

    // Load k-d tree
	KDTree* tree = new KDTree();
    if(!tree->load(filename)) mexErrMsgTxt("Error loading kd-tree\n");

    // return the program a pointer to the created tree
    plhs[0] = mxCreateDoubleScalar((long) tree);
}
#endif
