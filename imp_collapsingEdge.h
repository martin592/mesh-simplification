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
            a(pNull), b(pNull), collapsePoint(pNull), cost(-1) {}

        // Constructor
        collapsingEdge::collapsingEdge(const vector<UInt> & edge, costClass costObj)
        {
            UInt id_a = min(edge.at(0), edge.at(1));
            UInt id_b = max(edge.at(0), edge.at(1));

            a = grid->getNode(id_a);
            b = grid->getNode(id_b);

            pair<point,Real> optimumPair;

            // CostClass object for the cost computing
            optimumPair = costObj.getEdgeCost(edge);

            collapsePoint = optimumPair.first;
            cost = optimumPair.second;
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
            return (a.getId()==cE.getA() && b.getId()==cE.getB());
        }

        bool operator==(const vector<Uint> &edge)
        {
            return (a.getId()==edge.at(0) && b.getId()==edge.at(1));
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
                return (a.getId() < cE.getA().getId() );
        }

		private:

        void collapsingEdge::print()
        {
            cout<<"Edge with end points: "<< a.getId() <<"-"<< b.getId() <<" and cost "<< cost <<endl;
        }

}

#endif // IMP_COLLAPSINGEDGE_H_INCLUDED
