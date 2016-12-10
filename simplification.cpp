#include "simplification.h"

using namespace std;
using namespace geometry;

//
// Constructors
//
simplification<Triangle, MT, D>::simplification(mesh<Triangle,MT> *  _grid) :
    gridOperation(_grid), costObj(gridOperation), intersec(_grid), strucData(_grid)
{
    // Important control on coherence between the inputs:
    // error if the costType is GeoData but the meshType is GEO
    if (!(std::is_base_of<bcost<Triangle,MT,CostType>, CostType>::value))
                throw std::exception();

    // create the set of collapsing edges ordered by cost
    setupCollapsingSet();

    // define the fixed element
    findElemDontTouch();
}


void simplification<Triangle,MT,CostType>::setGrid(const mesh<Triangle,MT> * _grid)
{
    // call the setgrid from doctor
    gridOperation.setMeshPointer(_grid);
}


void simplification<Triangle,MT,CostType>::setupCollapsingSet()
{
    // variables
    vector<UInt> edge_i;
    pair<point, Real> costPair;
    vector<geoElement<Line>> edgeVector = gridOperation.getPointerToConnectivity()->getEdges();

    // loop on the edges list from the connections
    for (UInt i=0; i<edgeVector.size(); i++){

            // set the edge and the ids
            edge_i = edgeVector[i].getVertices();

            // take the cost information
            costPair = getCost(edge_i[0], edge_i[1]);
            ///costPair = getCost2(edge_i[0], edge_i[1]);  second technique to test, less controls

            // null cost means that there is no a valid collapse point
            if(costPair.second!=0.0)
            {
                // add the edge with the information cost to the collapsingSet and to collapseInfo
                collapsingSet.emplace(id1, id2, costPair.first, costPair.second);
                costObj.addCollapseInfo(id1, id2, costPair.first, costPair.second);
            }
    }
}



pair<point, Real> simplification<Triangle,MT,CostType>::getCost(const UInt id1, const UInt id2)
{
    // variables
    UInt id1, id2;
    Real c, c_min(0.0);
    point collapsePoint(pNull);
    vector<point> validPoints;

    id1 = edge[0];
    id2 = edge[1];

    validPoints = controlCollapse(edge);

    // loop on the valid points and keep the one with minimum cost
    if (validPoints.size()!=0){

                c_min = costObj.getCost(id1,id2,validPoints[0]);
                collapsePoint = validPoints[0];

                // loop on the valid points to take the one associated to the minimum cost
                for (UInt j=1; j<validPoints.size(); j++){
                    c = costObj.getCost(id1,id2,validPoints[j]);
                    if (c<c_min){
                        c_min = c;
                        collapsePoint = validPoints[j];
                    }
                }
    }

    // set the output pair with the collapse information
    pair<point, Real> costPair(collapsePoint,c_min);

    return(costPair);
}


pair<point, Real> simplification<Triangle,MT,CostType>::getCost2(const UInt id1, const UInt id2)
{
    // variables
    Real c;
    set<pair<point,Real>> costPairs;  /// TO ORDER WITH RESPECT TO THE SECOND ARGUMENT
    vector<point> pointsToTest;
    vector<UInt> edge={id1,id2};

    pointsToTest = costObj.getPointsList(id1,id2);
    costPairs.resize(pointsToTest.size());

    for (UInt j=0; j<pointsToTest.size(); j++){
            c=costObj.getCost(id1,id2,pointsToTest[j]);
            costPairs.insert(make_pair(pointsToTest[j],c));
    }

    for (auto it=costPairs.begin(); it!=costPairs.end(); it++){
        if(controlCollapsePoint(edge,it->first))  /// STILL TO IMPLEMENT
                return(*it);
    }
    return(make_pair(pNull,0));
}


void simplification<Triangle,MT,CostType>::findElemDontTouch()
{
      // variables
      point3d                             bar, barMesh(0.0,0.0,0.0);
      geoElementSize<Triangle>		      elem;
      set<geoElementSize<Triangle> >      barSet;

      // compute the barycenter
      for(UInt i=0; i<gridOperation.getPointerToMesh->getNumNodes(); ++i)
                            barMesh= barMesh + gridOperation.getPointerToMesh->getNode(i).getCoor();
      barMesh = barMesh/static_cast<Real>(gridOperation.getPointerToMesh->getNumNodes());

      // loop on the elements of the grid
      for(UInt i=0; i<gridOperation.getPointerToMesh->getNumElements(); ++i)
      {
        // set the element elem
        elem.setId(i);
        elem.setVertices(gridOperation.getPointerToMesh->getElement(i).getVertices());

        // find the barycenter of the specific element
        bar.reset();
        for (UInt j=0; j<3; j++)
            bar = bar + gridOperation.getPointerToMesh->getNode(elem.getVertices().at(j)).getCoor();
        bar = bar/3.0;

        // set the data size
        elem.setSize((bar-barMesh).norm2());

        // insert in the barSet
        barSet.insert(elem);
      }

      dontTouch =    true;
      elemDontTouchID = barSet.begin()->getId();
}


bool simplification<Triangle,MT,CostType>::controlLocalGrid(const vector<UInt> & involved, const vector<UInt> & edgePatch,
                                                  const vector<UInt> &toRemove)
{
      graphItem			      grafTmp;
      vector<UInt>            tmpEle;

      // create the graph of the involved elements
      graphItem grafInvolved(involved);

      // loop on the patch
      for(UInt i=0; i<edgePatch.size(); ++i)
      {
            grafTmp.clear();
            tmpEle.clear();

            // create the graph of the connections of node edgePatch[i]
            grafTmp.setConnected(connectivity.getNode2Elem(edgePatch[i]).getConnectedId);

            // remove the elements  which will be collapsed
            for(UInt s=0; s<toRemove.size(); ++s)
                        grafTmp.remove(toRemove[s]);

            // tmpEle contains the common elements between grafTmp and grafTmpNew
            set_intersection(grafTmp, grafInvolved, tmpEle);

            // if there are more or less common elements the control is failed
            if(tmpEle.size()!=2)
                        return(false);
            }

      return (true)
}

bool simplification<Triangle,MT,CostType>::controlDontTouch(const UInt id1, const UInt Id2)
{
    array<Uint>::iterator it1,it2;

    // take the vertices of the fixed element
    array<UInt,3> vertDontTouch = gridOperation.getPointerToMesh->getElem(elemDontTouchId).getVertices();

    it1=find(vertDontTouch.begin(),vertDontTouch.end(),id1);
    it2=find(vertDontTouch.begin(),vertDontTouch.end(),id2);

    // if I find one of the two end points in the array the control fails
    if(it1!=vertDontTouch.end() || it2!=vertDontTouch.end())
        return(false);

    // successful check
    return(true);
}

// specialize for geometric mesh
template<>
vector<point> simplification<Triangle,MT,CostType>::controlCollapse<MT::GEO>(const vector<UInt> & edge)
{
    //variables
    bool            valid=true;
    Real            area;
    UInt            id1, id2, oldP1;
    vector<UInt>    elemOnEdge, involved, allInvolved;
    vector<point3d> oldNormals;
    vector<point>   pointsToTest, validPoints;

    // set the ids
    id1 = edge[0];
    id2 = edge[1];

    // store the old point
    oldP1 = gridOperation.getPointerToMesh()->getNode(id1);
    // set the second id inactive
    gridOperation.getPointerToMesh()->setNodeInactive(id2);

    // set the involved mesh elements in the collapse
    elemsOnEdge = gridOperation.getElemsOnEdge();
    allInvolved = gridOperation.getElemsInvolvedInEdgeCollapsing(id1,id2);
    involved = gridOperation.getElemsModifiedInEdgeCollapsing(id1,id2);

    // set the normals before the test collapse
    oldNormals.resize(involved.size());
    for (UInt i=0; i<involved.size(); i++)
        oldNormals.emplace_back(gridOperation.getNormal(involved[i]));

    pointsToTest = costObj.getPointsList(id1,id2);
    validPoints.reserve(pointsToTest.size());

    auto conn = gridOperation.getPointerToConnectivity();
    auto oldConnections = conn->applyEdgeCollapse(id2, id1, elemsOnEdge, involved);

    for (UInt i=0; i<pointsToTest.size();i++)
    {
        valid=true;

        gridOperation.getPointerToMesh()->setNode(id1, pointsList[i]);

        strucData.update(involved);

		for (auto it = involved.cbegin(), auto oldNorm = oldNormals.cbegin(); it1 != involved.cend() && oldNorm != oldNormals.cend(); ++it, ++oldNorm)
        {
            // control inversions
            valid =  ((*oldNorm) * gridOperation.getNormal(*it) > sqrt(3)/2);
            if (!valid)
                break;

            // control on no self-intersections
            auto elems = strucData.getNeighbouringElements(*it);
            for (auto it2 = elems.cbegin(); it2 != elems.cend() && valid; ++it2)
					valid = !(intersec.intersect(*it, *it2));
            if (!valid)
                break;

            // control on no triangle degeneration
            area = gridOperation.getTriArea(*it);
            valid = !(area<TOLL);
            if (!valid)
                break;

            /// INSERT LOCAL CONTROL ON 2 ADJACENT EDGES?
        }

        if(valid)
            validPoints.push_back(pointsToTest[i]);
    }

    // restore connections
    conn->undoEdgeCollapse(id2, id1, oldConnections.first, oldConnections.second, elemsOnEdge, involved);

    // Restore list of nodes
    gridOperation.getPointerToMesh()->setNode(id1, oldP1);
    gridOperation.getPointerToMesh()->setNodeActive(id2);

    // Restore structured data
    strucData.update(involved);

    return(validPoints);
}

// specialize for mesh with data distrubuted oon
template<>
vector<point> simplification<Triangle,MT,CostType>::controlCollapse<MT::DATA>(const vector<UInt> & edge)
{
    //variables
    bool            valid=true;
    Real            area;
    UInt            id1, id2, oldP1;
    vector<UInt>    elemOnEdge, involved, allInvolved, toProject;
    vector<point3d> oldNormals;
    vector<point>   pointsToTest, validPoints;

    // set the ids
    id1 = edge[0];
    id2 = edge[1];

    // store the old point
    oldP1 = gridOperation.getPointerToMesh()->getNode(id1);
    // set the second id inactive
    gridOperation.getPointerToMesh()->setNodeInactive(id2);

    // store the elements of the mesh involved in the collapse
    elemsOnEdge = gridOperation.getElemsOnEdge();
    allInvolved = gridOperation.getElemsInvolvedInEdgeCollapsing(id1,id2);
    involved = gridOperation.getElemsModifiedInEdgeCollapsing(id1,id2);
    toProject = gridOperation.getDataModifiedInEdgeCollapsing(id1,id2);

    // fixed triangle control
    if(dontTouch)
            valid=controlDontTouch(edge[0], edge[1]);
    if(!valid)
        return(validPoints);

    // set the normals before the test collapse
    oldNormals.resize(involved.size());
    for (UInt i=0; i<involved.size(); i++)
        oldNormals.emplace_back(gridOperation.getNormal(involved[i]));

    // take the possible collapsePoints from the cost class object
    pointsToTest = costObj.getPointsList(id1,id2);
    validPoints.reserve(pointsToTest.size());

    auto conn = gridOperation.getPointerToConnectivity();
    auto oldConnections = conn->applyEdgeCollapse(id2, id1, elemsOnEdge, involved);

    // loop on the list of points to test
    for (UInt i=0; i<pointsToTest.size();i++)
    {
        valid=true;

        gridOperation.getPointerToMesh()->setNode(id1, pointsList[i]);

        strucData.update(involved);

        // Project data points and update data-element and element-data connections
        auto oldData = gridOperation.project(toProject, involved);

		for (auto it = involved.cbegin(), auto oldNorm = oldNormals.cbegin(); it1 != involved.cend() && oldNorm != oldNormals.cend(); ++it, ++oldNorm)
        {
            // control inversions
            valid =  ((*oldNorm) * gridOperation.getNormal(*it) > sqrt(3)/2 );
            if (!valid)
                break;

            // control on empty triangles
            valid = !(gridOperation.isEmpty(*it));
            if (!valid)
                break;

            // control on no self-intersections
            auto elems = strucData.getNeighbouringElements(*it);
            for (auto it2 = elems.cbegin(); it2 != elems.cend() && valid; ++it2)
					valid = !(intersec.intersect(*it, *it2));
            if (!valid)
                break;

            // control on no triangle degeneration
            area = gridOperation.getTriArea(*it);
            valid = !(area<TOLL);
            if (!valid)
                break;

            /// INSERT LOCAL CONTROL ON 2 ADJACENT EDGES?
        }

        if(valid)
            validPoints.push_back(pointsToTest[i]);

        gridOperation.undo(toProject,oldData);
    }

    // restore connections
    conn->undoEdgeCollapse(id2, id1, oldConnections.first, oldConnections.second, elemsOnEdge, involved);

    // Restore list of nodes
    gridOperation.getPointerToMesh()->setNode(id1, oldP1);
    gridOperation.getPointerToMesh()->setNodeActive(id2);

    // Restore structured data
    strucData.update(involved);

    return(validPoints);
}


void simplification<Triangle,MT,CostType>::refresh(){

      // variables
      UInt                              b, cont;
      point                             p;
      vector<UInt>	                    newId;
      array<Real,3>                     X;
      array<UInt,3>                     ids;
      geoElement<Triangle>              tria;
      vector<point>                     tmpPt;
      vector<geoElement<Triangle> >     tmpTr;

      // reserve memory space for the lists
      newId.reserve(gridOperation.getPointerToMesh->getNumNodes());
      tmpPt.reserve(gridOperation.getPointerToMesh->getNumNodes());
      tmpTr.reserve(gridOperation.getPointerToMesh->getNumElements());

      // loop on the nodes
      cont = 0;
      for(UInt i=0; i<gridOperation.getPointerToMesh->getNumNodes(); ++i)
      {
	    if(gridOperation.getPointerToMesh->getNode(i).active)
	    {
		    // take the points info from the mesh node
		    X = grid->getNode(i).getCoor();
		    b = gridOperation.getPointerToMesh->getNode(i).getBoundary();

		    // set the point p
		    p.setCoor(X);
		    p.setBoundary(b);
            p.setId(cont);

		    // new Id considers only the active nodes
		    newId[i] = cont;
		    ++cont;

		    // insert the point inside the temporary list
		    tmpPt.push_back(p);
	    }
      }

      // loop on the elements to create the temporary list of elements
      cont = 0;
      for(UInt i=0; i<gridOperation.getPointerToMesh->getNumElements(); ++i)
      {
	      // take the element info from the mesh
		  ids[0] = newId[gridOperation.getPointerToMesh->getElem(i).getVertices()[0]];
		  ids[1] = newId[gridOperation.getPointerToMesh->getElem(i).getVertices()[1]];
		  ids[2] = newId[gridOperation.getPointerToMesh->getElem(i).getVertices()[2]];
		  tmp = gridOperation.getPointerToMesh->getElem(i).getGeoId();

		  // create the new triangle
		  tria.setVertices(ids);
		  tria.setGeoId(tmp);
		  tria.setId(cont);
          ++cont;

		  // insert in the temporary list of elements
		  tmpTr.push_back(tria);
      }

      // clear and reset the lists
      gridOperation.getPointerToMesh->clear();
      for (UInt i=0; i<tmpPt.size())
        gridOperation.getPointerToMesh->insertNode(tmpPt.getCoor(), tmpPt.getBoundary());
      for (UInt i=0; i<tmpTr.size())
        gridOperation.getPointerToMesh->insertElem(tmpTr.getVertices(),tmpTr.getGeoId());

      // adjust the ids
      gridOperation.getPointerToMesh->setUpNodesIds();
      gridOperation.getPointerToMesh->setUpElemsIds();
}



void simplification<Triangle,MT,CostType>::update(const vector<UInt> & edge, const point collapsePoint)
{
    // variables
    Real                 c;
    UInt                 idCollapse = edge.at(0), id_ij;
    vector<UInt>   	     edgePatch, elemOnEdge;
    vector<vector<UInt>>  involvedEdges, newEdges;
    collapsingSet::iterator  it_collEdge;
    pair<point,Real>    costPair;

    involvedEdges.push_back(edge);

    // take the patch of the contracted edge (endpoints excluded)
    edgePatch = gridOperation.getNodesInvolvedInEdgeCollapsing(idCollapse, edge.at(1));
    elemOnEdge = gridOperation.getElemsOnEdge(idCollapse, edge.at(1));
    // loop on the nodes connected with the edge
    for (UInt i=0; i<edgePatch.size(); i++){

            // inner loop over all the nodes connected to the patch
            for (UInt j=0; j<gridOperation.getPointerToConnectivity->getNode2Node(edgePatch[i]).size();j++){

                id_ij = gridOperation.getPointerToConnectivity->getNode2Node(edgePatch[i].getId()).getConnected()[j];

                involvedEdges.emplace_back(array<UInt,2>({edgePatch[i].getId(), id_ij }));

                // the second endpoint is not anymore active
                if (id_ij==edge.at(1))
                    newEdges.emplace_back(array<UInt,2>({edgePatch[i].getId(), idCollapse }));
                else
                    newEdges.emplace_back(array<UInt,2>({edgePatch[i].getId(), id_ij }));
            }
    }

    // remove the involvedEdges from collapsingSet and cInfoList
    for (UInt i=0; i<involvedEdges.size();i++){

        // erase the edge from cInfoList and take the cost information
        c = costObj.eraseCollapseInfo(involvedEdges[i].at(0),involvedEdges[i].at(1));

        // search by cost in the collapsingSet
        it_collEdge = collapsingSet.find(c);

        // erase the edge in both the collapsing set
        collapsingSet.erase(it_colleEdge);
    }

    // update the connections in connectivity
    gridOperation.getPointerToConnectivity->applyEdgeCollapsing( edge.at(1), idCollapse, elemOnEdge, edgePatch);

    // update the cost class
    costObj.update(idCollapse);

    // adding the new edges to the collapsingSet and cInfoList computing their new costs
    for (UInt i=0; i<newEdges.size();i++){

                costPair = getCost(newEdges[i]);

                collapsingSet.emplace(id1, id2, costPair.first, costPair.second);
                costObj.addCollapseInfo(id1, id2, costPair.first, costPair.second);
    }
}

// specialize for geometric mesh
template<>
void simplification<Triangle,MT,CostType>::collapseEdge<MT::GEO>(const vector<UInt> & edge, const point collapsePoint){

    // variables
    UInt id1, id2;
    vector<UInt> toRemove, involved;

    id1 = edge.at(0).getId();
    id2 = edge.at(1).getId();

    // substitute the collapsePoint in the nodes list
    grid->getNode(id1).setCoor(collapsePoint);

    // set the flag active second end point to false
    grid->getNode(id2).setInactive();

    // remove the collapsed elements from the elems list in the mesh
    toRemove = getElemsOnEdge(id1,id2);
    for (UInt i=0; i<toRemove.size(); i++)
                    grid->eraseElem(toRemove[i]);

    // update the data structure
    involved = gridOperation.getElemsModifiedInEdgeCollapsing(id1,id2);
    strucData.update(involved);

}


// specialize for mesh with data distributed on
template<>
void simplification<Triangle,MT,CostType>::collapseEdge<MT::DATA>(const vector<UInt> & edge, const point collapsePoint){

    // variables
    UInt id1, id2;
    vector<UInt> toRemove, involved, oldData;

    id1 = edge.at(0).getId();
    id2 = edge.at(1).getId();

    // substitute the collapsePoint in the nodes list
    grid->getNode(id1).setCoor(collapsePoint);

    // set the flag active second end point to false
    grid->getNode(id2).setInactive();

    // remove the collapsed elements from the elems list in the mesh
    toRemove = getElemsOnEdge(id1,id2);
    for (UInt i=0; i<toRemove.size(); i++)
                    grid->eraseElem(toRemove[i]);

    // update the data structure
    involved = gridOperation.getElemsModifiedInEdgeCollapsing(id1,id2);
    strucData.update(involved);

    // project the data
    toProject = gridOperation.getDataModifiedInEdgeCollapsing(id1,id2);
    auto oldData = gridOperation.project(toProject, involved);
}


void simplification<Triangle,MT,CostType>::simplificate(const UInt numNodesMax)
{
      // variables
      UInt              counter=0;
      UInt	            numNode=grid->getNumNodes();
      UInt              numNodeStart=grid->getNumNodes();
      time_t 	        start, end;
      Real              dif;
      collapsingEdge::iterator    minCostEdge;
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
      while(numNode>numNodesMax)
      {
            // take the first valid collapsing edge with the minimum cost
            minCostEdge = collaspingSet.begin();

            p = minCostEdge->getCollapsePoint();
            edge[0] = minCostEdge->getA();
            edge[1] = minCostEdge->getB();

            // collapse the edge in collapsePoint p
            collapseEdge<MT>(edge, p);

            // updates collapsingSet, collapseInfo and the connections
            update(edge, collapsePoint);

            // decrease the number of nodes
            --numNode;
      }

      // reset and updates the lists of the mesh
      refresh();

      time(&end);
      dif = difftime(end,start);
      cout << "The mesh size passed from " << numNodeStart << " to " << grid->getNumNodes() << " nodes" << endl;
      cout << "Simplification process completated in " <<  dif << " s" << endl;
}
