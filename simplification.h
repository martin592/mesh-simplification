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
#include "doctor.hpp"
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
            - The collapse function must be enriched to take in account the case the edge
              is close to the border
     */

    template<typename SHAPE, typename MT> class simplification : public doctor<SHAPE,MT>
    {
    };


    /*! Specialization for the triangular meshes */

    template<> class simplification<Triangle,MT> : public doctor<Triangle,MT>
{
      //
      // Attributes
      //

      public:

		  /*! Set of edges ordered by cost (ascending order) */
		  set<collapsingEdge>      collapsingSet;


      private:

          /*! Variable responsible for the cost computing of the edges */
          costClass        costObj;

    //
    // Constructors
    //

    public:
		  /*! (Default) Constructor
                \param _meshPointer pointer to the mesh */
		  simplification(smart_ptr<mesh<Triangle,MT>> _grid = nullptr);

		  /*! Method which changes the pointer to the mesh
		    \param _meshPointer pointer to the mesh */
		  void setGrid(const smart_ptr<mesh<Triangle,MT>> _grid);

          /*! Method that builds the set of CollapsingEdge ordered by cost.
              The method uses the edge list from the connections and adds the cost information. */
		  void setUpCollapsingSet();

		  /*! Method which makes the refresh of connections and other variables after collapse */
          void refresh();

	//
	// Method which creates the list
	//

	public:

		  /*! Method which updates the collapsingSet and the connections after each contraction */
		  void update(const vector<UInt> & edge, const point collapsePoint, const vector<UInt> & involved);

	//
	// Methods for the controls
	//

	public:
		 /*! Control on the correctness of the collapse.
		 The method involves:
              - check of the inversion of the normals
              - control of the maintaimnent of two triangles insisting on each edge
		     \param edge contracted
		     \param pNew new point */
		 bool control(const vector<UInt> & edge, const point pNew);

	//
	// Methods which make the simplification
	//

	public:
	    /*! Modification of the mesh after edge contraction
            //  In particular it performs:
            //     -  the update of the new node in the list with the id of the first end point,
            //     -  the inactivation of the second end point in the nodes list. */
        void simplification<SHAPE,MT>::collEdge(const vector<UInt> & edge, const point collapsePoint);   /// DA METTERE IN DOCTOR?

        /*! Method which iteratively contracts the edge with minimum cost until reaching a maximum
              amount of nodes
		      \param numNodesMax maximum number of nodes */
        void simplificate(const UInt numNodesMax);

}

#endif // SIMPLIFICATION_H_INCLUDED
