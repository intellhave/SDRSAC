#include "mex.h"
#include "KDTree.h"

void KDTree::from_matlab_matrix(const mxArray *matstruct){
    /// Retrieves datapoints
    {
        /// Retrieve memory
        mxArray* points_mex = mxGetField(matstruct, 0, "points");
        double* points_data = mxGetPr(points_mex);        
        
        /// Retrieve dataset size
        this->ndim = mxGetN(points_mex);
        this->npoints = mxGetM(points_mex);
        this->points.resize(npoints, Point(ndim));
        // cout << "ndim" << ndim << endl;
        // cout << "npoints" << npoints << endl;
        
        /// Fill memory
        for( int i=0; i<npoints; i++ )
            for( int j=0; j<ndim; j++ )
                points[i][j] = points_data[ i + j*npoints ];    
    }
    
    /// Retrieve nodes
    {
        /// Retrieve memory
        mxArray* nodes_mex = mxGetField(matstruct, 0, "nodes");
        double* nodes_data = mxGetPr(nodes_mex);
        
        /// Retrieve data size
        int datasize = mxGetM(nodes_mex);
        // mexPrintf("size(nodes_mex) = %d\n",datasize);
        nodesPtrs.resize( datasize, NULL );
        
        /// Fill memory
        for(int i=0,off=0; i<nodesPtrs.size(); i++,off+=4){
            assert(i<nodesPtrs.size());
            assert(off+4<datasize);
            
            nodesPtrs[i] = new Node();            
            nodesPtrs[i]->LIdx = nodes_data[off+0];
            nodesPtrs[i]->RIdx = nodes_data[off+1];
            nodesPtrs[i]->pIdx = nodes_data[off+2];
            nodesPtrs[i]->key  = nodes_data[off+3];
        }
    }
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){
    if(nrhs!=1) mexErrMsgTxt("varargin{1}  is the struct pointer");
    if(nlhs!=1) mexErrMsgTxt("varargout{1} is the tree pointer");
    
    /// Instantiate tree
    KDTree* tree = new KDTree();
    tree->from_matlab_matrix(prhs[0]);
   
    /// Store pointer in matlab
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    double* pointer_to_tree = mxGetPr(plhs[0]);
    pointer_to_tree[0] = (intptr_t) tree;
}
