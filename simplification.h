#ifndef SIMPLIFICATION_H_INCLUDED
#define SIMPLIFICATION_H_INCLUDED


#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <functional>
#include <exception>
#include <type_traits>
#include <numeric>
#include <iterator>
#include <utility>

#include "utility.hpp"
#include "bmeshOperation.h"
#include "bcost.hpp"
#include "collapsingEdge.h"
#include "intersection.hpp"
#include "structuredData.hpp"
#include "geoPoint.hpp"
#include "geoElement.hpp"
#include "geoElementSize.h" /// STILL TO DO

namespace geometry
{

    /*! Class for the simplification process on the mesh elements by iterative edge collapse.
        The selection of the edge to contract is done with respect to a cost functional
        which measures the loss of information for the analysis. In case of a meshType MT::GEO,
        the loss of information is just geometrical, otherwise for meshType MT::DATA the cost
        functional takes into account the loss of statistical information as well.

        At each iteration the edge with minimum cost is collapsed into a predefined point.
        The edges are organized in a set in ascending order with respect to their costs called
        collapsingSet which is used to extract the edge to contract. In case the edge has to be
        searched for its ending points, refer to the unordered set cInfoList contained in the
        cost class object.

        The point for the collapse is chosen among the two ending points of the selected edge,
        their middle point and the so-called optimal point. The latter derives from an explicit
        analytical formula based on proposition 5.2.2 from the thesis essay "Advanced Techniques for the
        Generation and the Adaptation of Complex Surface Meshes" by F. Dassi.

        The selection of the collapsing point depends on several controls. First of all, the control
        in the cost class takes into account wheter the edge or an endpoint belongs to the border.
        In addition, the elected point has to pass the following geometric controls:
            - the collapse does not invert triangles
            - the contraction has to not bring to mesh self-intersections
            - no degenerated triangles are generated
        In case the meshType is of type DATA, there are the further data controls.
            - the fixed element is not affected by the edge collapse
            - no empty triangle are created after the data projection

        Each edge collapse is implemented by setting the id of the collapse point to the id of an
        end point (the one with minor id) and setting to inactive the flag of the second end point.
        The elements which insist on the edge are removed from the mesh.
        Subsequently, the costs of the involved edges are recomputed and the connections are
        locally updated. In case of mesh DATA the involved data points are projected.

        The class presents three templates:
            1. SHAPE, geometrical shape of the elements of the mesh
            2. MeshType, GEO or DATA, whether the mesh contains only mesh points or mesh and data points
            3. CostType, OnlyGeo or GeoData, wheter the cost used is based only on the geometrical information
                                             or also on the statical
        The critical situation corresponds to a CostType::GeoData in combination with a MeshType::GEO. This
        conflicting situation is handled at compile time in the constructor checking the templates coherence
        with type_trades.

        The class contains the following attributes:
            \param gridOperation, object of conditional class depending on the meshType:
                        - if MT::GEO, it corresponds to a meshInfo object
                        - if MT::DATA, it is a projection object to handle distrubuted data
            \param costObj, cost class of costType CT
            \param collapsingSet, set of collapsingEdges ordered by cost in ascending order
            \param intersec, interesection object for the related control
            \param structData, data structure necessary to support the intersection control
            \param dontTouch, boolean to indicate if the fixed element is used
            \param elemDontTouchId, id of the fixed element
     */

    template<typename SHAPE, MeshType MT, CostType CT> class simplification
    {
    };


    /*! Specialization for the triangular meshes */

    template<>
    class simplification<Triangle, MT, CostType>
{
      //
      // Attributes
      //
      private:

        /*! MeshOperation Object  */
        bmeshOperation<Triangle,MT>      gridOperation;

        /*! CostClass object for the cost computing of the edges */
        costType                  costObj;

        /*! Set of edges ordered by cost (ascending order) */
        set<collapsingEdge>      collapsingSet;

        /*!  object for the control of triangle intersections */
        intersection<Triangle>    intersec;

        /*! object for the bounding boxes structure  */
        structuredData<Triangle>     strucData;

        /*! fixed element to not touch */
        //  boolean to indicate if it used or not
        bool          dontTouch;
        // id of the element
        UInt          elemDontTouchId;

    //
    // Constructors
    //

    public:

		  /*! (Default) Constructor
                \param _grid pointer to the mesh */
		  simplification(mesh<Triangle,MT> * _grid = nullptr);

		  /*! Method which changes the pointer to the mesh
		    \param _grid pointer to the mesh */
		  void setGrid(const mesh<Triangle,MT> * _grid);

          /*! Method that builds the set of CollapsingEdge ordered by cost and
              the unordered set of cInfoList in the cost class
              The method uses the edge list from the connections and adds the cost information. */
		  void setupCollapsingSet();

		  /*!   IMPORTANT
              Method which takes the cost data for the contraction of the edge.
              The function operates by step:
                - extracts from the cost class object the list of possible collapse points
                - controls the validity of the points
                - takes from the cost class object the minimum cost value
                \param id1, id2, ids of the end points of teh edge to contract
                The method returns a pair with the minimum cost and the collapse point associated*/
		  pair<point,Real> getCost(const UInt id1, const UInt id2);

		  /*!   IMPORTANT
                Method that differently to the previous version first compute all the costs for
                the possible collapse points and then check if the one with minimum cost is valid.
                If not the point with the second smaller cost is checked and so on.
                This implementation should be more efficient because saves time for the controls*/
		  pair<point,Real> getCost2(const UInt id1, const UInt id2);


		  /*! Method which makes the refresh of connections and other variables after collapse */
          void refresh();

    //
    // Methods for the elemDontTouch handling
    //

        /*! get/set for dontTouch*/
        INLINE void activeDontTouch(){
		      dontTouch=true;
		      };

        INLINE void disactiveDontTouch(){
		      dontTouch=false;
		      };

        INLINE bool getDontTouch(){
		      return(dontTouch);
		      };

        /*! get/set for elemDontTouch */
        INLINE void setElemDontTouch(UInt _elemDontTouch){
            elemDontTouch=_elemDontTouch;
            };

        INLINE UInt getElemDontTouch(){
            return(elemDontTouch);
            };


       /*! Method which automatically finds the element to fix */
		  void findElemDontTouch();


      //
      // Methods for specific controls
      //

      private:

        /*! Method which controls that each edge maintains exactly two adjacent triangles */
        bool controlLocalGrid(const vector<UInt> & involved, const vector<UInt> & edgePatch,
                                                                        const vector<UInt> & toRemove);

        /*! Method to check that the selected edge do not edge affect the element to not touch
		      \param toRemove, elements which insist on the collapsing edge
		      \param involved, elements which change after the contraction
            It return true if the fixed element is not touched, false otherwise*/
        bool controlDontTouch(const vector<UInt> & toRemove, const vector<UInt> & involved);

    //
    // Global Controls
    //
    private:
        /*! Method for check on the validity on the collapse of the edge into the specific input point P
            It is a template method on the meshType.
            If MT::DATA, there are the additional controls on the empty triangles and on the fixed element;
            in this case, the method projects the data points involved and undoes the projection after each control.
		    \param edge to contract
		    \param P collapse point */

        template<MT>
        bool controlCollapsePoint(const vector<UInt> & edge, const point P)

        /*! Method to check the validity on the collapse of the edge. It considers all the possible collapse
            points proposed by the cost class.
            The function tries the update of connections, performs the controls and undo the operations.
            It is a template method on the meshType.
            If MT::DATA, there are the additional controls on the empty triangles and on the fixed element;
            in this case, the method projects the data points involved and undoes the projection after each control.
		    \param edge to contract */
        template<MT>
		bool controlCollapse(const vector<UInt> & edge);

	//
	// Method fort the update of the mesh and the connections
	//

	private:

		  /*! Method which updates the mesh, the connectivities, the data structure, the collapsingSet
		  and the cInfoList after each contraction.
		  It is a template method on the meshType.
		  In case mesh::DATA, there is the further update of the distribution of the data points.
		  \param edge
		  \param collapsePoint */
		  template<MT>
		  void update(const vector<UInt> & edge, const point collapsePoint);

	//
	// Methods which make the simplification
	//

	public:

        /*! Method which iteratively contracts the edge with minimum cost until reaching a maximum
              amount of nodes
		      \param numNodesMax maximum number of nodes */
        void simplificate(const UInt numNodesMax);


//        /*! Method which iteratively contracts the edge with minimum cost until reaching a maximum
//              amount of nodes
//		      \param numNodesMax maximum number of nodes
//		      \param collapseNumber number of edge of contractions per each iteration */
//        void simplificateGreedy(const UInt numNodesMax, const UInt collapseNumber);

}

#endif // SIMPLIFICATION_H_INCLUDED
