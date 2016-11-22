#ifndef COSTCLASS_H_INCLUDED
#define COSTCLASS_H_INCLUDED

#include <vector>
#include <cassert>
#include "point.hpp"

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

class costClass
{

    public:

		  /*! Vector of the geometric matrices Q for each node */
		  vector<vec>           Qvec;

          /*! Connections */
		  connect<Triangle,MT>   connectivity;

      //
      // Constructors
      //

      public:
          /*! Constructor */
          costClass(connect<Triangle,MT>   _conn);

          /*! Copy constructor */
          costClass(costClass & cC);

          /*! Destructor */
          ~costClass();

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

          /*! Update of the Q matrix of nodeInvolved and node nodeCollapse
                \param  nodeIds, list of nodes  */
          void updateQVector(const vector<UInt> nodeIds);

      //
      // Methods for finding matrices and point for the collapse and cost computing
      //

      public:
          /*! Method that returns the optimal point from Qtilde
                \param edge */
          point createOptimalPoint(const vector<UInt> & edge) const;

		  /*! Method which creates the list of nodes to test
              The method considers inversion problems, the optimal point exceptions and the border end points
                \param edge
                \param newNodes list of test points */
		  void createPointList(const vector<UInt> & edge, const vector<point> & newNodes) const;

		  /*! Method which returns the collapse cost of an edge
                \param QEdge geometrical matrix for the edge
                \param pNew collapse point */
		  Real getEdgeCost(const vec & QEdge, const point pNew) const;

		  /*! Method which return the collapse point with minimum cost and the cost itself
                \param edge */
 		  pair<point, Real> getEdgeCost(const vector<UInt> & edge) const;

}

/*! Include implementations of class members. */
#include "implementation/imp_costClass.h"


#endif // COSTCLASS_H_INCLUDED
