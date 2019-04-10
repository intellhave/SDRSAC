/// OLD SOURCE FOR CLOSEST POINT

/// Computes the closest point in the set to the query point "p"
int KDTree::closest_point(Point p){
    double neigh_dst = 0;
    return closest_point(p, neigh_dst);
}
int KDTree::closest_point(Point p, double& neigh_dst){
    // search closest leaf
    int dim = 0; // variable to cycle through dimensions
    Node* leaf = nodesPtrs[0];
    for (;;) {
        // Is leaf node... this is my stop
        if( leaf->pIdx >= 0 )
            break;

        // Not a leaf... browse through
        if( p[dim] <= leaf->key )
            leaf = nodesPtrs[ leaf -> LIdx ];
        else
            leaf = nodesPtrs[ leaf -> RIdx ];
        dim = (dim + 1) % ndim;
    }

    double cdistsq = distance_squared( p, points[leaf->pIdx] ); 	// best distance at the moment
    // cout << "first approximation: " << leaf->idx+1 << endl;
    int closest_neighbor = leaf->pIdx;
    check_border_distance(ROOT, 0, p, cdistsq, closest_neighbor); 		//check if anything else can do better
    
    neigh_dst = sqrt( cdistsq );
    return closest_neighbor;
}

/** @see closest_point
 *
 * This function is the algorithm in support of "closest_point" for
 * closest point computation.
 *
 * @param nodeIdx the index of the node to check for the current recursion
 * @param dim     the dimension to check in the current recursion
 * @param pnt     the query point
 * @param cdistsq the euclidean distance for query to the point "idx"
 * @param idx     the index to the "currently" valid closest point
 */
void KDTree::check_border_distance(int nodeIdx, int dim, Point pnt, double& cdistsq, int& idx){
    Node* node = nodesPtrs[ nodeIdx ];
    // cout << "checking node: " << node->idx+1 << endl;

    // Are we at a leaf node? check if condition and close recursion
    if( node->pIdx >= 0 ){
        // is the leaf closer in distance?
        float dsq = distance_squared(pnt, points[ node->pIdx ] );
        if (dsq < cdistsq){
            cdistsq = dsq;
            idx = node->pIdx;
            // cout << "updated optimal with: " << node -> idx+1 << endl;
        }
        return;
    }

    // The distance squared along the CURRENT DIMENSION between the point and the key
    float ndistsq = (node->key - pnt[dim])*(node->key - pnt[dim]);
    // cout << "distance to key: " << ndistsq << " optimal current distance: " << cdistsq << "(dim:"<<dim<<")"<<endl;

    // If the distance squared from the key to the current value is greater than the
    // nearest distance, we need only look in one direction.
    if (ndistsq > cdistsq) {
        if (node->key > pnt[dim])
            check_border_distance(node->LIdx, (dim+1)%ndim, pnt, cdistsq, idx);
        else
            check_border_distance(node->RIdx, (dim+1)%ndim, pnt, cdistsq, idx);
    }
    // If the distance from the key to the current value is less than the nearest distance,
    // we still need to look in both directions.
    else {
        //cout << "both directions need to be checked" << endl;
        check_border_distance(node->LIdx, (dim+1)%ndim, pnt, cdistsq, idx);
        check_border_distance(node->RIdx, (dim+1)%ndim, pnt, cdistsq, idx);
    }
}
