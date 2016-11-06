/*! \file dataPoint.hpp
	\brief Declaration of a class representing a point with an associated datum. */
	
#ifndef HH_DATAPOINT_HH
#define HH_DATAPOINT_HH

#include "point.hpp"

namespace geometry
{
	using namespace std;
	
	/*! Class inheriting point and expanding it adding the datum associated to the point. */
	class dataPoint final : public point
	{
		private:
			/*! The datum. */
			Real datum;
			
		public:
			//
			// Constructors
			//
			
			/*! Constructor. 
				\param c	array with coordinates 
				\param ID	point id
				\param dat	datum */
			dataPoint(const array<Real,3> & c, const UInt & ID = 0, const UInt & dat = 0.);
			
			/*! Constructor.
				\param p	point
				\param dat	datum */
			dataPoint(const point & p, Real const & dat);
			
			/*! Synthetic copy constructor.
				\param p	another point */
			dataPoint(const dataPoint & p) = default;
			
			//
			// Operators
			//
			
			/*! Copy assignment operator.
				\param p	another point */
			dataPoint & operator=(const dataPoint & p);
			
			//
			// Get methods
			//
			
			/*! Get the datum.
				\return 	the datum */
			inline Real getDatum() const {return datum;};
			
			//
			// Set methods
			//
			
			/*! Set the datum.
				\param dat	the new datum */
			inline void setDatum(Real const & dat) {datum = dat;};
			
		private:
			/*! Print to output the point data.
				\param out	the output string */
			virtual void print(ostream & out) const;
	}; 
}

#endif
