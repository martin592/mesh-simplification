#ifndef INLINE_COLLAPSINGEDGE_H_INCLUDED
#define INLINE_COLLAPSINGEDGE_H_INCLUDED

namespace geometry{

    //
    // Get methods
    //
    public:
        /*! Get the first end point */
        INLINE UInt getA() const {
            return a;
        };

        /*! Get the second end point */
        INLINE UInt getB() const {
            return b;
        };

        /*! Get the collapse point */
        INLINE point getCollapsePoint() const {
            return collapsePoint;
        };

        /*! Get the cost */
        INLINE Real getCost() const {
            return cost;
        };

    //
    // Set methods
    //
    public:
        /*! Set the first end point*/
        INLINE void setA(const UInt a_) {
            a = a_;
        };

        /*! Set the second end point*/
        INLINE void setB(const UInt b_) {
            b = b_;
        };

        /*! Set the optimal point*/
        INLINE void setCollapsePoint(const point cp) {
            collapsePoint = cp;
        };

        /*! Set the cost*/
        INLINE void setCost(const Real c) {
            cost = c;
        };
}

#endif // INLINE_COLLAPSINGEDGE_H_INCLUDED
