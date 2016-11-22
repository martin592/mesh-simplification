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

costClass::costClass(connect<Triangle,MT>  _conn)
{
    connectivity = _conn;

    setupQVector();
}


costClass::costClass(costClass & cC){
    Qvec = cC.QVec;
    connectivity = cC.connectivity;
}

~costClass(){
    Qvec.clear();
    connectivity.clear();
}

//
// Preparatory methods for the geometric matrices
//

vec costClass::createK_p(const UInt nodeId, const UInt elemId)
    {

    assert(nodeId<grid->getNumNodes());
    assert(elemId<grid->getNumElements());

    // variables
    Real    noteTerm;
    point   normal;
    vec     K_p;

    //	take the normal
    normal = getNormal(elemId);

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


vec costClass::createQ(const UInt nodeId){

        assert(nodeId< grid->getNumNodes());

        // variables
        vec	    QTmp, K_p;

        // initialize the 4x4 matrix to 0
        QTmp.assign(16,0.0);

        // take the connected nodes to nodeId
        for(UInt i=0; i<connectivity.Node2Elem[nodeId].size()); ++i)
        {
            // compute the K_p matrix
            K_p = createK_p(nodeId, connectivity.node2elem[nodeId].getConnected().at(i));

            // update the matrix Q
            QTmp = QTmp + K_p;

        }

        // return the matrix
        return(QTmp);

}


vec costClass::createQ(const vector<UInt> & edge)
    {
        // initialize the 4x4 matrix to 0
        vec QEdge;
        QEdge.assign(16,0.0);

        // sum of the Q matrices of the end points
        QEdge = Qvec[edge.at(0)]+Qvec[edge.at(1)];

        return(QEdge)
}


void costClass::setupQVector()
    {
    Qvec.clear();
    Qvec.resize(grid->getNumNodes());

    // loop on the nodes of the mesh to create all the matrices
    for(UInt i=0; i<grid->getNumNodes(); ++i)
        Qvec[i] = createQ(i);
}


void costClass::updateQVector(const vector<UInt> nodeIds)
{
    for (UInt i=0; i<nodeIds.size(); i++)
            Qvec[nodeIds[i]] = createQ(nodeIds[i]);
}

//
// Cost computing
//

point costClass::createOptimalPoint(const vector<UInt> & edge)
{
      // variables
      UInt           cont=0;
      point			 pNew;
      vec	         QTmp;

      Epetra_SerialDenseMatrix  QTilde;
      Epetra_SerialDenseVector	   v,b;
      Epetra_SerialDenseSolver	solver;

      // construct the edge matrix Qtmp
      QTmp = createQTilde(edge);

      // remove the elements I do not care for the solver
      QTmp[12] = 0.0;
      QTmp[13] = 0.0;
      QTmp[14] = 0.0;
      QTmp[15] = 1.0;

      // create Qtilde
      QTilde.Shape(4,4);
      for(UInt i=0; i<4; ++i)
            for(UInt j=0; j<4; ++j, ++cont)
                QTilde(i,j) = QTmp[cont];

      v.Size(4);
      b.Size(4);

      // set the note term b
      b(0) = 0.0;
      b(1) = 0.0;
      b(2) = 0.0;
      b(3) = 1.0;

      // set the variables of the solver
      solver.SetMatrix(QTilde);
      solver.SetVectors(v,b);

      double res = solver.Solve();

      // no solution so I return pNull
      if(res==-1.0)
                return(pNull);

      // otherwise v contains the coordinates of the optimal point
      pNew.setX(v(0));
      pNew.setY(v(1));
      pNew.setZ(v(2));

      return(pNew);
}


void costClass::createPointList(const vector<UInt> & edge,const vector<point> & newNodes)
{
   // variables
    point     	pMid, pOpt;
    point p1    =	grid->getNode(edge.at(0));
    point p2 	=	grid->getNode(edge.at(1));
    UInt bound1 = grid->getNode(edge.at(0)).getBoundary();
    UInt bound2 = grid->getNode(edge.at(1)).getBoundary();
    vector<UInt>	geoIds1, geoIds2;
    vector<point>	newNodeTmp;

    // reserve on the nodes-to-test list
    newNodeTmp.clear();
    newNodeTmp.reserve(5);

    // switch on different border conditions
    switch(bound1)
    {
      case(0):
	      switch(bound2)
	      {
		// both intern points
		case(0):

			// compute the optimal point
			pOpt = createPointFromMatrix(edge);
			pOpt.setBoundary(0);

			if(pOpt!=pNull)
                newNodeTmp.push_back(pOpt);

			// compute the middle point
			pMid.replace(p1,p2,0.5);
			pMid.setBoundary(0);

			newNodeTmp.push_back(pMid);

			newNodeTmp.push_back(p1);
			newNodeTmp.push_back(p2);

			break;
		// p2 on the border
		case(1):
			newNodeTmp.push_back(p2);
			break;
	      }
	      break;
      case(1):
	      switch(bound2)
	      {
		// p1 on the border
		case(0):
			newNodeTmp.push_back(p1);
			break;
		// p1 e p2 both on the bordo
		/// ma si deve controllare davvero se sono effettivamente un edge (??)
		case(1):
			if(isBoundary(edge))
			{
			    // control on the geoIds around to preserve the mesh angles
			    geoIdAround(edge.at(0), &geoIds1);
			    geoIdAround(edge.at(1), &geoIds2);

			    if((geoIds1.size()>2) && (geoIds2.size()==2))	              /// GEOIDS SIZE??
                    newNodeTmp.push_back(p1);
			    else if((geoIds1.size()==2) && (geoIds2.size()>2))
                        newNodeTmp.push_back(p2);
                     else
                    {
                        // compute the middle point
                        pMid.replace(p1,p2,0.5);
                        // this time it is on th border
                        pMid.setBoundary(1);
                        newNodeTmp.push_back(pMid);

                        newNodeTmp.push_back(p1);
                        newNodeTmp.push_back(p2);
			    }
			}
			break;
	      }
	      break;
    }

    // control on the identified nodes to check if they may lead to the inversion problem
    newNodes.clear();
    newNodes.reserve(newNodeTmp.size());
    for(UInt i=0; i<newNodeTmp.size(); ++i)
            if(controlColl(edge, newNodeTmp[i]))                   ///CONTROLCOLL?
                newNodes.push_back(newNodeTmp[i]);
}


Real costClass::getEdgeCost(const vec & QEdge, const point pNew){
      // variables
      vector<Real>	 v, vTmp;

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
      return(inner_product(v1.begin(), v1.end(), vTmp.begin(), 0.0));
}


pair<point, Real> costClass::getEdgeCost(vector<UInt> & edge){

      // variables
      Real 	  	            cost, costTmp;
      point                 pNew;
      vector<point>	        newNodes;
      vec	                QEdge;
      pair<point, Real >    result;

      // create the edge matrix QEdge
      QEdge = createQ(edge);

      // create the list of test nodes
      createPointList(edge, &newNodes);
l
      // no valid test nodes
      if(newNodes.size()==0)
      {
            result.first = pNull;
            return(result);
      }

      // take the first pair point-cost and then compare with others in newNodes
      result.first = newNodes[0];
      result.second = getEdgeCost(QEdge, newNodes[0]);

      for(UInt i=1; i<newNodes.size(); ++i)
      {
            // compute the cost
            costTmp = getEdgeCost(QEdge, newNodes[i]);

            // substitute if the cost is minor
            if(costTmp<result.second)
            {
                result.first = newNodes[i];
                result.second = costTmp;
            }
      }

      // return the pair <collapsePoint, cost>
      return(result);
}


#endif // IMP_COSTCLASS_HPP_INCLUDED
