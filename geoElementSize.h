#ifndef GEOELEMENTSIZE_H_INCLUDED
#define GEOELEMENTSIZE_H_INCLUDED

#include <cassert>
#include <iostream>

#include "geoElement.hpp"

namespace geometry
{

using namespace std;

/*! Class which inherits from geoElement and adds a certain data "size" to the mesh element.
GeoElementSize will be used to sort elements in base of the value size and employ that for the
simplification algorithm.*/

template<typename SHAPE> class geoElementSize : public geoElement<SHAPE>
{
	public:

		/*! Additional value associated to the element */
		Real size;
	//
	// Constructors
	//
	public:
		/*! Costruttore */
		geoElementSize();

		/*! Costruttore di copia */
		geoElementSize<SHAPE> (const geoElementSize<SHAPE> &E);

		/*! Operatore uguaglianza */
		geoElementSize<SHAPE> operator=(const geoElementSize<SHAPE> &E);
	//
	// Operatori
	//
	public:
		/*! Operatore non uguaglianza */
		bool operator!=(const geoElementSize<SHAPE> &E) const;

		/*! Operatore non uguaglianza */
		bool operator==(const geoElementSize<SHAPE> &E) const;

		/*! Operatore minore */
		bool operator<(const geoElementSize<SHAPE> &E) const;

	//
	// Set/get
	//
	public:
		/*! Funzione che permette di accedere al size dell'elemento */
		inline Real getGeoSize() const;

		/*! Funzione che permette di settare l'id geometrico */
		inline void setGeoSize(const Real & _size);

	//
	// Stampa
	//
	public:
		/*! Print to screen */
		void print();
};

//-------------------------------------------------------------------------------------------------------
// IMPLEMENTATION
//-------------------------------------------------------------------------------------------------------

//
// Constructors
//
template<typename SHAPE> geoElementSize<SHAPE>::geoElementSize() : geoElement<SHAPE>::geoElement()
{
	size = 0.0;
}

template<typename SHAPE> geoElementSize<SHAPE>::geoElementSize(const geoElementSize<SHAPE> &E) : geoElement<SHAPE>(E)
{
	size = E.size;
}

template<typename SHAPE> geoElementSize<SHAPE> geoElementSize<SHAPE>::operator=(const geoElementSize<SHAPE> &E)
{
	this->geoElement<SHAPE>::operator=(E);
	size  = E.size;
	return *this;
}

//
// Get/Set
//
template<typename SHAPE> inline Real geoElementSize<SHAPE>::getGeoSize() const
{
	return(size);
}

template<typename SHAPE> inline void geoElementSize<SHAPE>::setGeoSize(const Real & _size)
{
	size= _size;
}

//
// Operators
//

template<typename SHAPE> bool geoElementSize<SHAPE>::operator!=(const geoElementSize<SHAPE> &E) const
{
	if(E.size==this->size) 	return(this->graphItem::operator!=(E));
	else			return(true);

}


template<typename SHAPE> bool geoElementSize<SHAPE>::operator==(const geoElementSize<SHAPE> &E) const
{
	return(this->geoElement<SHAPE>::operator==(E));

}


template<typename SHAPE> bool geoElementSize<SHAPE>::operator<(const geoElementSize<SHAPE> &E) const
{
	if(this->size==E.size) 	return(this->graphItem::operator<(E));
	else			return(this->size<E.size);
}



template<typename SHAPE> void geoElementSize<SHAPE>::print()
{
	cout << "Id         : " << this->getId() << endl;
	cout << "GeoId      : " << this->getGeoId() << endl;
	cout << "Size       : " << size << endl;
	cout << "Num Points : " << this->getNumVertices() << endl;
	cout << "Nodes Id's : ";

	for(UInt i=0; i<this->getNumVertices(); ++i)
	{
		cout << geoElement<SHAPE>::getVertices.at(i) << " ";
	}
	cout << endl;
}

}

#endif // GEOELEMENTSIZE_H_INCLUDED
