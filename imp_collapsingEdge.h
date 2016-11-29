#ifndef IMP_COLLAPSINGEDGE_H_INCLUDED
#define IMP_COLLAPSINGEDGE_H_INCLUDED

namespace geometry
{
	using namespace std;

	//
	// Constructors
	//

        // Default Constructor
        collapsingEdge::collapsingEdge() :
            a(-1), b(-1), collapsePoint(pNull), cost(-1) {}

        // Constructor
        collapsingEdge::collapsingEdge(const vector<UInt> & edge, costClass costObj)
        {
            a = min(edge.at(0), edge.at(1));
            b = max(edge.at(0), edge.at(1));

            pair<point,Real> optimumPair;

            // CostClass object for the cost computing
            optimumPair = costObj.getEdgeCost(edge);

            collapsePoint = optimumPair.first;
            cost = optimumPair.second;
        }

        collapsingEdge::collapsingEdge(const UInt a_, const UInt b_, const point collPt, const Real c)
        {
            a = a_;
            b = b_;

            collapsePoint = collPt;
            cost = c;
        }


	//
	// Operators overload
	//

        void operator=(const collapsingEdge &cE)
        {
            a = cE.a;
            b = cE.b;
            collapsePoint = cE.collapsePoint;
            cost = cE.cost;

        }


        bool operator==(const collapsingEdge &cE)
        {
            return (a==cE.getA() && b==cE.getB());
        }

        bool operator==(const vector<Uint> &edge)
        {
            return (a==edge.at(0) && b==edge.at(1));
        }

        bool operator!=(const collapsingEdge &cE)
        {
            return (!(*this==cE));
        }

        bool operator!=(const vector<Uint> &edge)
        {
            return (!(*this==edge));
        }

        bool operator<(const collapsingEdge &cE)
        {
            if (cost!=cE.cost)
                return (cost < cE.getCost());
            else
                return (a < cE.getA() );
        }

		private:

        void collapsingEdge::print()
        {
            cout<<"Edge with end points: "<< a <<"-"<< b <<" and cost "<< cost <<endl;
        }

}

#endif // IMP_COLLAPSINGEDGE_H_INCLUDED
