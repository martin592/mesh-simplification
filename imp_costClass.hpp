#ifndef IMP_COSTCLASS_HPP_INCLUDED
#define IMP_COSTCLASS_HPP_INCLUDED

// operator overload for the sum between vector<double>
vec operator+(const vec & v1, const vec & v2)
{
    assert(v1.size()==v2.size());
    vec v;
    v.resize(v1.size());

    for (auto i=0; i<v1.size(); i++){
        v[i] = v1[i] + v2[i];
    }
    return(v);
}

// operator overload for the product between vector<double> column*row = matrix
vec operator*(const vec & v1, const vec & v2)
{
    vec v;
    for (auto i=0; i<v1.size(); i++){
        for (auto j=0; j<v2.size(); j++)
            v.push_back(v1[i]*v2[j]);
    }
    return(v);
};


//
// Constructor
//

costClass<Meshtype>::costClass(meshInfo<Triangle,MT>   _gridInfo)
{
    gridInfo= _gridInfo;

    setupQVector();
}


costClass<Meshtype>::costClass(costClass & cC){
    Qvec = cC.QVec;
    gridInfo = cC.gridInfo;
}

~costClass(){
    Qvec.clear();
    gridInfo.connectivity.clear();
}


vector<vec> costClass<Meshtype>::getQVec()
{
        return(Qvec);
}

void costClass<Meshtype>::setQVec(vector<vec> & _Qvec)
{
    Qvec = _Qvec;
}

void costClass<Meshtype>::updateQVector(const vector<UInt> nodeIds)
{
    for (UInt i=0; i<nodeIds.size(); i++)
            Qvec.at(nodeIds[i]) = createQ(nodeIds[i]);
}

void costClass<Meshtype>::updateCostObject(const connect<MeshType> _conn, const vector<UInt> nodeIds)
{

    gridInfo.connectivity = _conn;

    updateQVector( nodeIds);
}


void costClass<Meshtype>::addQ(vec & _Q)
{
    Qvec.push_back(_Q);
}

void costClass<Meshtype>::getQ(UInt id)
{
    return (Qvec.at(id))
}


//
// Preparatory methods for the geometric matrices
//

vec costClass<Meshtype>::createK_p(const UInt nodeId, const UInt elemId)
    {

    assert(nodeId<grid->getNumNodes());
    assert(elemId<grid->getNumElements());

    // variables
    Real    noteTerm;
    point3D   normal;
    vec     K_p;

    //	take the normal
    normal = gridInfogetNormal(elemId);

    // set the noteTerm
    noteTerm = (-1.0)*(grid->getNode(nodeId)*normal);

    // create an extended vector of the normal
    vec normal_ext={normal.getX(), normal.getY(), normal.getZ(), noteTerm};

    // initialize the matrix to 0
    K_p.assign(16,0.0);

    K_p = normal_ext*normal_ext;

    // return the matrix
    return(K_p);

}


vec costClass<Meshtype>::createQ(const UInt nodeId){

        assert(nodeId< grid->getNumNodes());

        // variables
        vec	    QTmp, K_p;

        // initialize the 4x4 matrix to 0
        QTmp.assign(16,0.0);

        // take the connected nodes to nodeId
        for(UInt i=0; i<gridInfo.connectivity.getNode2Elem(nodeId).size()); ++i)
        {
            // compute the K_p matrix
            K_p = createK_p(nodeId, gridInfo.connectivity.getNode2Elem(nodeId).getConnected().at(i));

            // update the matrix Q
            QTmp = QTmp + K_p;

        }

        // return the matrix
        return(QTmp);

}


vec costClass<Meshtype>::createQ(const vector<UInt> & edge)
    {
        // initialize the 4x4 matrix to 0
        vec QEdge;
        QEdge.assign(16,0.0);

        // sum of the Q matrices of the end points
        QEdge = Qvec[edge.at(0)] + Qvec[edge.at(1)];

        return(QEdge)
}


void costClass<Meshtype>::setupQVector()
    {
        Qvec.clear();
        Qvec.resize(grid->getNumNodes());

        // loop on the nodes of the mesh to create all the matrices
        for(UInt i=0; i<grid->getNumNodes(); ++i)
            Qvec[i] = createQ(i);

}


Real costClass<Meshtype>::getCost(const vec & QEdge, const point pNew){

      // variables
      vec	 v, vTmp;
      Real      c;

      // v contains the coordinates of pNew and 1 in fourth position
      v.assign(4,1.0);
      for(UInt i=0; i<3; ++i)
                v[i] = pNew[i];

      // vTmp is the result of multiplication between QEdge and v
      vTmp.assign(4, 0.0);
      vTmp[0] = QEdge.at(0)*v[0] + QEdge.at(4)*v[1] + QEdge.at(8)*v[2]  + QEdge.at(12)*v[3];
      vTmp[1] = QEdge.at(1)*v[0] + QEdge.at(5)*v[1] + QEdge.at(9)*v[2]  + QEdge.at(13)*v[3];
      vTmp[2] = QEdge.at(2)*v[0] + QEdge.at(6)*v[1] + QEdge.at(10)*v[2] + QEdge.at(14)*v[3];
      vTmp[3] = QEdge.at(3)*v[0] + QEdge.at(7)*v[1] + QEdge.at(11)*v[2] + QEdge.at(15)*v[3];

      // return value is the inner product
      c = inner_product(v.begin(), v.end(), vTmp.begin(), 0.0);

      return(c);
}


#endif // IMP_COSTCLASS_HPP_INCLUDED
