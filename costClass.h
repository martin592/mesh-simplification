#ifndef COSTCLASS_H_INCLUDED
#define COSTCLASS_H_INCLUDED

#include <vector>
#include <cassert>
#include "point.hpp"
#include "meshInfo.hpp"

/// include for the optimal point computation

#include "Epetra_ConfigDefs.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_LinearProblem.h"


using vec = vector<Real>;

/*! The class provides the value of the cost functional associated to an edge.
    It contains the following:
    <ol>
    <li> the vector of geometric matrices Q (4x4) for each node, reshaped as vector (16x1)
    </ol>
    Currently the cost functional is purely geometric.
    The implementation for the complete cost functional has to be done.
*/

template<Meshtype MT = MeshType::GEO>
class costClass<MT>
{

    protected:

		  /*! Vector of the geometric matrices Q for each node */
		  vector<vec>           Qvec;

          /*! Connections */
		  meshInfo<Triangle,MT>   gridInfo;

      //
      // Constructors
      //

      public:

          /*! Constructor */
          costClass(meshInfo<Triangle,MT>   _gridInfo);

          /*! Copy constructor */
          costClass(costClass & cC);

          /*! Destructor */
          ~costClass();

      //
      // Methods get/set/update
      //
          /*! Method which returns the vector of Q matrices */
		  vector<vec> getQVec() const;

          /*! Method which provides the vector of Q matrices to the class
                \param  _Qvec, input vector of Qs  */
		  void setQVec(vector<vec> & _Qvec);

          /*! Update of the entire costClass object
                \param  _conn, connect object
                \param  nodeIds, list of nodes to update */
		  void costClass<Meshtype>::updateCostObject(const connect<MeshType> _conn, const vector<UInt> nodeIds)


          /*! Method which adds a Q matrix to the vector QVec
                \param  _Q, matrix to add to the list Qvec */
		  void addQ(vec & _Q);

          /*! Method which take the Q matrix corresponding to nodeId i
                \param  id, node id  */
		  void getQ(UInt id);

      //
      // Methods for the update
      //

          /*! Update of the Q matrix of nodeInvolved and node nodeCollapse
                \param  nodeIds, list of nodes to update
              Pay attention: this method requires that the connectivity node2elem is updated in gridInfo.*/
          void updateQVector(const vector<UInt> nodeIds);

          /*! Update of the entire costClass object, including the connectivity of gridInfo
                \param  _conn, connect object
                \param  nodeIds, list of nodes to update */
		  void updateCostObject(const connect<MeshType> _conn, const vector<UInt> nodeIds)

      //
      // Methods for the building of the geometric matrices associated to nodes and edges of the mesh
      //

      public:

		  /*! Method which builds matrix K_p
                \param nodeId node identifier
                \param elemId element identifier */
		  vec createK_p(const UInt nodeId, const UInt elemId) const;

		  /*! Method which returns the (reshaped) matrix Q for a specific node
		        \param nodeId  node identifier */
		  vec createQ(const UInt nodeId);

          /*! Method which builds the geometric matrix Q for the edge
                \param edge pointer to the edge */
		  vec createQ(const vector<UInt> & edge) const;

		  /*! Method which sets the vector of Q matrices of the nodes */
		  void setupQVector();

      //
      // Method for cost computing
      //

		  /*! Method which returns the collapse cost of an edge
                \param QEdge geometrical matrix for the edge
                \param pNew collapse point */
		  Real getCost(const vec & QEdge, const point pNew) const;

}

/*! Include implementations of class members. */
#include "implementation/imp_costClass.h"


#endif // COSTCLASS_H_INCLUDED
