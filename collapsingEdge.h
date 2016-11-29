#ifndef COLLAPSINGEDGE_H_INCLUDED
#define COLLAPSINGEDGE_H_INCLUDED

#include <iostream>
#include <utility>

#include "geoElement.hpp"
#include "point.h"
#include "costClass.h"


namespace geometry{

using namespace std;

/*! Edge class thought specifically for an iterative simplification process which collapses the edge with minimum cost.
    It contains the following:
        <ol>
        <li> ids of the end points;
        <li> optimal point for the collapse;
        <li> cost value.
        </ol>
    The contraction information derives from the costClass methods, given as input in the collapsingEdge constructor.
*/

class collapsingEdge
{

    private:

    // Ids of the two ending points of the edge
    UInt a;
    Uint b;

    // Point for the collapse
    point collapsePoint;

    // Edge cost
    Real cost;

	//
	// Constructors
	//

    public:

        /*! Default Constructor */
        collapsingEdge();

        /*! Constructor from simple edge
                \param edge, without cost information
                \param costObj costClass obkject */
        collapsingEdge(const vector<UInt> & edge, costClass costObj);


        collapsingEdge(const UInt a_, const UInt b_, const point collPt, const Real c);

        /*! Copy constructor
            \param collapsing edge to copy  */
        collapsingEdge(const collapsingEdge & cE) = default;

        /*! Destructor */
        ~collapsingEdge() = default;

	//
	// Operators overload
	//

    public:
        /*! Assigning operator */
        & collapsingEdge operator=(const collapsingEdge &cE);

        /*! Identity operator */
        bool operator==(const collapsingEdge &cE) const;

        /*! Identity operator for edge input */
        bool operator==(const vector<UInt> &edge) const;

        /*! Inequality operator */
        bool operator!=(const collapsingEdge &cE)

        /*! Inequality operator for edge input */
        bool operator!=(const vector<UInt> &edge) const;

        /*! Minor operator */
        bool operator<(const collapsingEdge &cE) const;

    //
    // Get methods
    //

    public:
        /*! Get the id of the first end point */
        point getA() const;

        /*! Get the id of the second end point */
        point getB() const;

        /*! Get the collapse point */
        point getCollapsePoint() const;

        /*! Get the cost */
        point getCost() const;

    //
    // Set methods
    //

    public:
        /*! Set the id of first end point*/
        void setA(const point a_);

        /*! Set the id of second end point*/
        void setB(const point b_);

        /*! Set the optimal point*/
        void setCollapsePoint(const point cp);

        /*! Set the cost*/
        void setCost(const point c);

    //
    // Print methods
    //

    public:
        /*! Print to output the  edge with the cost*/
        void print() const;

};

}

#ifdef INLINED
#include "inline/inline_collapsingEdge.hpp"
#endif

/*! Include implementations of class members. */
#include "implementation/imp_collapsingEdge.h"


#endif // COLLAPSINGEDGE_H_INCLUDED
