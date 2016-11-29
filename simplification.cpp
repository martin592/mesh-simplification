#include "simplification.h"

using namespace std;
using namespace geometry;

//
// Constructors
//
simplification<Triangle,MT>::simplification() : meshOperation<Triangle,MT>()
{
}

simplification<Triangle,MT>::simplification(mesh<Triangle,MT> *  _grid) : meshOperation<Triangle,MT>(_grid)
{
    costObj->setgridInfo = gridOperation.connectivity;

    // create the set of collapsing edges ordered by cost
    setupCollapsingSet();
}


void simplification<Triangle,MT>::setGrid(const mesh<Triangle,MT> * _grid)
{
    // call the setgrid from doctor
    meshOperation.setMeshPointer(_grid);

    // create the set of collapsing edges ordered by cost
    setupCollapsingSet();
}


void simplification<Triangle,MT>::setupCollapsingSet()
{
    vector<UInt> edge_i;
    pair<point, Real> costPair;

    costObj->setupQVector();

    vector<geoElement<Line>> edgeVector = gridOperation.connectivity.getEdges();

    // loop on the edges list from the connections
    for (UInt i=0; i<edgeVector.size(); i++){

            edge_i[0] = edgeVector[i].at(0);
            edge_i[1] = edgeVector[i].at(1);

            costPair = getEdgeCost(edge_i);

            collapsingSet.emplace(edge_i.at(0), edge_i.at(1), costPair.first, costPair.second);

            gridOperation.connectivity.edgeSet.at(i).setGeoSize(costPair.second);  /// STILL TO IMPLEMENT AS vector<GeoElementSize<Line>>
    }
}


//
// Cost computing
//

point simplification<Triangle,MT>::createOptimalPoint(const vector<UInt> & edge)
{

      // variables
      UInt           cont=0;
      point			 pNew;
      vec	         QTmp;

      Epetra_SerialDenseMatrix  QTilde;
      Epetra_SerialDenseVector	   v,b;
      Epetra_SerialDenseSolver	solver;

      // construct the edge matrix Qtmp
      QTmp = createQ(edge);

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


void simplification<Triangle,MT>::createPointList(const vector<UInt> & edge, const vector<point> & newNodes)
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
			pOpt = createOptimalPoint(edge);
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

			    if((geoIds1.size()>2) && (geoIds2.size()==2))
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
            if(controlCollapse(edge, newNodeTmp[i]))
                newNodes.push_back(newNodeTmp[i]);
}


pair<point, Real> simplification<Triangle,MT>::getEdgeCost(vector<UInt> & edge){

      // variables
      Real 	  	            cost, costTmp;
      point                 pNew;
      vector<point>	        newNodes;
      vec	                QEdge;
      pair<point, Real >    result;

      // create the edge matrix QEdge
      QEdge = costObj->createQ(edge);

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
      result.second = costObj->getCost(QEdge, newNodes[0]);

      for(UInt i=1; i<newNodes.size(); ++i)
      {
            // compute the cost
            costTmp = costObj->getCost(QEdge, newNodes[i]);

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


void simplification<Triangle,MT>::refresh(){

      // variables
      UInt                              tmp, cont, id1, id2, id3;
      Real                              x, y, z;
      point                             p;
      geoElement<Triangle>              tria;
      vector<point>                     tmpPt;
      vector<geoElement<Triangle> >     tmpTr;
      vector<UInt>	                    newId;


      // reserve memory space for the lists
      newId.reserve(grid->getNumNodes());
      tmpPt.reserve(grid->getNumNodes());
      tmpTr.reserve(grid->getNumElements());

      // set to true doing a loop over the elements
      for(UInt i=0; i<grid->getNumElements(); ++i)
      {
	    if(!isTriangleDegenerate(i))
	    {
		  id1 = grid->getElem(i).getConnectedId(0);
		  id2 = grid->getElem(i).getConnectedId(1);
		  id3 = grid->getElem(i).getConnectedId(2);

		  grid->getNode(id1).setActive();
          grid->getNode(id2).setActive();
          grid->getNode(id3).setActive();
	    }
      }

      // loop on the nodes
      cont = 0;
      for(UInt i=0; i<grid->getNumNodes(); ++i)
      {
	    if(grid->getNode(i).active)
	    {
		    // take the points info from the mesh node
		    x   = grid->getNode(i).getX();
		    y   = grid->getNode(i).getY();
		    z   = grid->getNode(i).getZ();
		    tmp = grid->getNode(i).getBoundary();

		    // set the point p
		    p.setX(x);	  p.setY(y);	p.setZ(z);
		    p.setBoundary(tmp);         p.setId(cont);

		    // new Id considers only the active nodes
		    newId[i] = cont;
		    ++cont;

		    // insert the point inside the temporary list
		    tmpPt.push_back(p);
	    }
      }

      // loop on the elements to create the temporary list of elements
      cont = 0;
      for(UInt i=0; i<grid->getNumElements(); ++i)
      {
	    if(!isTriangleDegenerate(i))
	    {
	      // take the element info from the mesh
		  id1 = newId[grid->getElem(i).getConnectedId(0)];
		  id2 = newId[grid->getElem(i).getConnectedId(1)];
		  id3 = newId[grid->getElem(i).getConnectedId(2)];
		  tmp = grid->getElem(i).getGeoId();

		  // create the new triangle
		  tria.setConnectedId(0,id1);
		  tria.setConnectedId(1,id2);
		  tria.setConnectedId(2,id3);
		  tria.setGeoId(tmp);
		  tria.setId(cont);
          ++cont;

		  // insert in the temporary list of elements
		  tmpTr.push_back(tria);
	    }
      }

      // clear and reset the lists
      grid->clear();
      grid->insertNode(&tmpPt);
      grid->insertElement(&tmpTr);

      // adjust the ids
      grid->setupIds();

      // setUp
      setup();                   /// CHE SETUP E'?? DELLE CONNESSIONI?

}



void simplification<Triangle,MT>::update(const vector<UInt> & edge, const point collapsePoint, const vector<UInt> & involved)
{
    // variables
    UInt                 idCollapse = edge.at(0), id_ij;
    vector<UInt>   	     edgePatch;
    vector<vector<UInt>>  involvedEdges, newEdges;
    collapsingSet::iterator  it;

    involvedEdges.push_back(edge);

    // remove the collapsed elements from the elems list in the mesh
    toRemove = getElemsOnEdge(id1,id2);
    for (UInt i=0; i<elems.size(); i++)
                    grid->eraseElem(toRemove[i]);

    // take the patch of the contracted edge (endpoints excluded)
    edgePatch = getNodesInvolvedInEdgeCollapsing(idCollapse, edge.at(1));

    // loop on the nodes connected with the edge
    for (UInt i=0; i<edgePatch.size(); i++){

            // inner loop over all the nodes connected to the patch
            for (UInt j=0; j<gridOperation.connectivity.getNode2Node(edgePatch[i]).size();j++){

                id_ij = gridOperation.connectivity.getNode2Node(edgePatch[i]).getConnected().at(j);

                involvedEdges.emplace_back(array<UInt,2>({edgePatch[i], id_ij }));

                // the second endpoint is not anymore active
                if (id_ij==edge.at(1))
                    newEdges.emplace_back(array<UInt,2>({edgePatch[i], idCollapse }));
                else
                    newEdges.emplace_back(array<UInt,2>({edgePatch[i], id_ij }));
            }
    }

    // remove the involvedEdges from the collapsingSet
    for (UInt i=0; i<involvedEdges.size();i++){
        it = collapsingSet.find(involvedEdges[i]);
        if(it!=collapsingSet.end())
                collapsingSet.erase(it);
    }

    // update the connections in connectivity
    gridOperation.connectivity.applyEdgeCollapsing( edge.at(1), idCollapse, toRemove, edgePatch);

    // update the vector of Q matrices with the one related to the new point and update the connectivity in costObj
    costObj->updateCostObject(gridOperation.connectivity, edgePatch);

    // adding the new edges to the collapsingSet computing their costs
    for (UInt i=0; i<newEdges.size();i++)
            collapsingSet.emplace(newEdges[i]);
}



void simplification<Triangle,MT>::simplificate(const UInt numNodesMax)
{
      // variables
      UInt              counter=0;
      UInt	            numNode=grid->getNumNodes();
      UInt              numNodeStart=grid->getNumNodes();
      time_t 	        start, end;
      Real              dif;
      collapsingEdge    minCostEdge;
      point 	        p;
      vector<UInt>		edge, elemToUpdate;

      // control if I am already below
      if(numNodesMax>=grid->getNumNodes())
      {
            cout << "The number of mesh points is " << grid->getNumNodes();
            cout << ", already below to the given threshold " << numNodesMax << endl;
            return;
      }

      // print
      cout << "Simplification process..." << endl;
      time(&start);

      // iterative collapse until numNodeMax is reached
      while(numNode>numNodesMax && counter<numNodeStart)
      {
            // take the first collapsing edge with the minimum cost
            minCostEdge = collapsingSet[0];

            p = minCostEdge.getCollapsePoint();
            edge[0] = minCostEdge.getA();
            edge[1] = minCostEdge.getB();

            if(control(edge, p))
            {
                // collapse the edge in collapsePoint p
                collapseEdge(&edge, p);

                // updates collapsingSet and the connections
                update(edge, collapsePoint, involved);

                // decrease the number of nodes
                --numNode;
            }
            else
            {
                collapsingSet.erase(0);
            }

            // increase the counter variable for each simplification trial (even if the control is failed!)
            ++counter;

            cout << numNode << " nodes with a maximum allowed of " << numNodesMax << "        \r";
      }

      // make a final refresh updating the lists from the mesh
      refresh();

      time(&end);
      dif = difftime(end,start);
      cout << "The mesh size passed from " << numNodeStart << " to " << grid->getNumNodes() << " nodes" << endl;
      cout << "Simplification process completated in " <<  dif << " s" << endl;
}
