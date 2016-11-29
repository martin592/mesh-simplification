#ifndef SIMPLIFICATION_H_INCLUDED
#define SIMPLIFICATION_H_INCLUDED


#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <functional>
#include <numeric>

#include "utility.hpp"
#include "meshOperation.h"
#include "costClass.hpp"
#include "collapsingEdge.h"



using vec = vector<Real>;

namespace geometry
{

    /*! Class for the simplification process on the mesh elements by iterative edge collapse.
        The selection of the edge to contract is done with respect to a cost functional
        which measures the loss of geometrical (and statistical) information for the analysis;
        at each iteration the edge with minimum cost is collapsed into a predefined point.
        The edges are organized in a set in ascending order with respect to their costs.

        The point for the collapse is chosen among the two ending points of the selected edge,
        their middle point and the so-called optimal point. The latter derives from an explicit
        formula based on proposition 5.2.2 from the thesis essay "Advanced Techniques for the
        Generation and the Adaptation of Complex Surface Meshes" by F. Dassi.
        The selection of the collapsing point takes into account wheter the edge or a point of its
        lies on the border.

        Each edge collapse is implemented by setting the id of the collapse point to the id of an
        end point (the one with minor id) and setting to inactive the flag of the second end point.
        The elements which insist on the edge are removed from the mesh.
        Subsequently, the costs of the involved edges are recomputed and the connections are
        locally updated.

        TO DO:
            - The statistical part of the simplication process has to be implemented yet
     */

    template<typename SHAPE, typename MT> class simplification
    {
    };


    /*! Specialization for the triangular meshes */

    template<> class simplification<Triangle,MT>
{
      //
      // Attributes
      //

      public:

		  /*! Set of edges ordered by cost (ascending order) */
		  set<collapsingEdge>      collapsingSet;

      private:

          /*! MeshOperation Object  */
          meshOperation<Triangle,MT>      gridOperation;

          /*! Pointer to a costClass object for the cost computing of the edges */
          costClass  *                  costObj;

    //
    // Constructors
    //

    public:

		  /*! (Default) Constructor
                \param _meshPointer pointer to the mesh */
		  simplification(mesh<Triangle,MT> * _grid = nullptr);

		  /*! Method which changes the pointer to the mesh
		    \param _meshPointer pointer to the mesh */
		  void setGrid(const mesh<Triangle,MT> * _grid);

          /*! Method that builds the set of CollapsingEdge ordered by cost.
              The method uses the edge list from the connections and adds the cost information. */
		  void setupCollapsingSet();

		  /*! Method which makes the refresh of connections and other variables after collapse */
          void refresh();


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

          /*! Method which return the collapse point with minimum cost and the cost itself
                \param edge */
 		  pair<point, Real> getEdgeCost(const vector<UInt> & edge) const;


	//
	// Method fort the update of the mesh and the connections
	//

	public:

		  /*! Method which updates the collapsingSet and the connections after each contraction */
		  void update(const vector<UInt> & edge, const point collapsePoint, const vector<UInt> & involved);


	//
	// Methods which make the simplification
	//

	public:

        /*! Method which iteratively contracts the edge with minimum cost until reaching a maximum
              amount of nodes
		      \param numNodesMax maximum number of nodes */
        void simplificate(const UInt numNodesMax);

}

#endif // SIMPLIFICATION_H_INCLUDED
