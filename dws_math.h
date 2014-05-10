/**************************************************************************/
// FILE NAME: dwsmath.h
// VERSION: $Revision: 1.3 $ 
// REMARKS: Dylan Schwilk's math template library
//------------------------------------------------------------------------*/
// Copyright © 2002 - 2005 */
/*------------------------------------------------------------------------*/
// AUTHOR/S:  Dylan W. Schwilk */
//------------------------------------------------------------------------*/
// CREATION DATE: 020830
//--------+------+-------------------------------------+------------------*/
// REVISION LOG:
// $Log: dws_math.h $
// Revision 1.3  2005/02/18 04:49:30  dylan
// Fixed deleted char.
//
// Revision 1.2  2005/01/26 17:31:51  schwilkd
// modified code to meet ansi standards and compile on mingw, cygwin,
// linux. Added Doxygen comments and removed non-templated math functions
// to dws_numerical library (dws_numerical.h and dws_numerical.cpp).  The
// dws_math library now consists only of the template header file
// dws_math.h.  I also un-capitalized the utility functions (abs, etc) to
// conform with my coding standards.
//
// Revision 1.1  2005/01/26 17:26:21  schwilkd
// Initial revision
//
//                   
/******************************************************************************/

//////////////////////////////////////////////////////////////////////////////
// GNU
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Library General Public License as published 
// by the Free Software Foundation; either version 2 of the License, or 
// (at your option) any later version.  This library is distributed in the
// hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//////////////////////////////////////////////////////////////////////////////


/*!
 * \file
 * The dws_math libary provides a set of templated mathematical objects and 
 * functions. 
 */ 

#ifndef _DWS_MATH_H
    #define _DWS_MATH_H
    
    #include <cmath>
    #include <limits>
//    #include <iterator>
    
    namespace DWS
    {
    
    /*! \defgroup dws_math DWS math template library
     * @{
     */
     
        // conventional utility functions -- templated
        // These three all should be available in std::, but I keep running into
        // problems, hence, my own versions.
        template <class T>
        inline T abs(T t){return ( t > 0 ? t : -t);}
    
        template <class T>
        inline T min(T a, T b){return (a < b ? a : b);}

        template <class T>
        inline T max(T a, T b){ return (a > b ? a : b);}
        
        template <class T>
        inline T median(T a, T b, T c){
            return ( max(a,max(b,c))==a ? max(b,c) : (c>b ? max(a,b):max(a,c) ));
        }
    
        template <class T>
        inline T sign(T a){return (a > 0 ? 1 : -1);}
    
        template <class T>
        inline T sqr(T a){return (a*a);} //!< square
        
        
////////////////////////////////////////////////////////////////////////////////   
//       Structs and classes
////////////////////////////////////////////////////////////////////////////////
    
        // Struct: TPos
        // ---------
        /*! Position type
         * 
         * Templated 2D position type useful for graphics and grids.  
         * Instaniate with either floating or integral types.
         */
        template <class T> class TPos
        {
            public:
            TPos() : _x(0), _y(0) {}
            TPos(T x, T y) : _x(x),_y(y) {}
            T _x, _y;
    
            bool operator<(const TPos<T>& other) const;  //!< Tests less than 
            bool operator==(const TPos<T>& other) const; //!< Tests equality
            T operator-(const TPos<T>& other) const; //!< Subtract two positions: gives euclidean distance
        };
    
        // TPos functions
        template <class T>
        bool TPos<T>::operator< (const TPos<T>& other) const
        {return(_x < other._x || (!(other._x < _x) && _y < other._y));}
    
        template <class T>
        bool TPos<T>::operator==(const TPos<T>& other) const
        {return (_x==other._x && _y==other._y);}
    
        template <class T>
        T TPos<T>::operator-(const TPos<T>& other) const
        {return ( std::sqrt(sqr(_x-other._x) + sqr(_y-other._y)) );}
    
    
    
        // class TLine
        // ----------
        /*! Templated line class useful for graphics.  Instaniate with either 
         *  floating or integral types.
         */
        template <class T> class TLine
        {
            public:
            TLine() : _s(0,0),_f(0,0) {}
            TLine(TPos<T> s, TPos<T> f) : _s(s), _f(f) {}
            TPos<T> _s, _f;
            bool operator<(const TLine<T>& other) const;  /// test less_than 
            bool operator==(const TLine<T>& other) const; /// Test equality of two lines
        };
    
    
        // TLine functions
        template <class T>
        bool TLine<T>::operator<(const TLine<T>& other) const
        {return(_s<other._s || (!(other._s<_s) && _f<other._f));}
    
        template <class T>
        bool TLine<T>::operator==(const TLine<T>& other) const
        {return (_s==other._s && _f==other._f);}
    
        
        /* Type RangeT
         * ----------- */
         /*! Type for storing a high-low range.
          * 
          * This templated class is for use with any floating point or integer
          * type.  A RangeT object supports a variety of operations such as
          * sets and unions.  A RangeT contains two public data members
          * \arg _min and \arg _max.
          */
        template <class T> class TRange
        {
            public:
            T _min;
            T _max;
            TRange() {_min = 1;_max = -1;}
            TRange(T x, T y) {_min = DWS::min(x,y);_max = DWS::max(x,y);}
            template <class A> TRange(const TRange<A>& copy):_min(copy._min),_max(copy._max){}
            bool Empty() const {return (_min > _max);}
            void Add(T n)
            {
                if(Empty()) _min = _max = n;
                else
                {
                    _min = std::min(_min,n);
                    _max = std::max(_max,n);
                }
            }
            void MakeEmpty() {_min = 1;_max = -1;}
            bool Singleton() const {return (_min == _max);}
        };  // TRange
    
    
        // Function GetRange(R1, R2)
        // -------------------------
        /*! Returns range of the medians of input ranges.
         * (In the MacClade manual, this is the circle with the X function).
         */
        template <class T>
        TRange<T> GetRange(const TRange<T>& r1, const TRange<T>& r2)
        {
            TRange<T> result;
            double y,z; // use double arithmatic internally
                        // FIXME: is this necessary?
    
            y = median(r1._min, r1._max, r2._min);
            z = median(r2._max, r2._min, r1._max);
            result._min = std::min(y,z);
            result._max = std::max(y,z);
            return result;
        }
    
        // Function: Intersect(r1,r2)
        // --------------------------
        /*! Returns the intersection of two ranges.
         */
        template <class T>
        TRange<T> Intersect(const TRange<T>& r1,const TRange<T>& r2)
        {
            TRange<T> result;
            if(r1.Empty() || r2.Empty())
            {
                result._min = 1;
                result._max = -1;
            }
            else
            {
                result._min = DWS::max(r1._min, r2._min);
                result._max = DWS::min(r1._max, r2._max);
            }
            return result;
        }
    
        // Function: Union
        // ---------------
        /*! Union of two ranges.
         */
        template <class T>
        TRange<T> Union(const TRange<T>& r1,const TRange<T>& r2)
        {
            TRange<T> result;
            result._min = DWS::min(r1._min, r2._min);
            result._max = DWS::max(r1._max, r2._max);
            return result;
        }
    
        // Average
        // -------
        /*! Mean of a range.
         */
        template <class T>
        TRange<T> RangeAverage(const TRange<T>& r)
        {
            TRange<T> result;
            result._min = result._max = (r._max+r._min) / 2.0;
            return result;
        }
        
    /*!
      @}
      */
        
    }; /* namespace DWS */
    
#endif
// _DWS_MATH_H
