/*************************************************************************/
// FILE NAME:dws_numerical.h
// VERSION: $Revision: 1.2 $
// FUNCTION: Dylan Schwilk's numerical function library
// REMARKS:
//------------------------------------------------------------------------
// Copyright © 2003, 2004, 2005
//------------------------------------------------------------------------
// AUTHOR/S:  Dylan W. Schwilk 
//------------------------------------------------------------------------
// CREATION DATE: 19/04/03 23:56
//--------+------+---------------------------------------------+----------
// REVISION LOG: 
//
// $Log: dws_numerical.h $
//
/*************************************************************************/

//////////////////////////////////////////////////////////////////////////
// GNU This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public License
// as published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.  This library
// is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General
// Public License for more details.  You should have received a copy
// of the GNU Library General Public License along with this library;
// if not, write to the Free Software Foundation, 59 Temple Place -
// Suite 330, Boston, MA 02111-1307, USA.
//////////////////////////////////////////////////////////////////////////

/*! \file The dws_numerical libary provides a set of non-templated
 * mathematical functions.  Many of these functions are modeled after
 * those in <em> Numerical Recipes in C</em>.
 */ 

#ifndef _DWS_NUMERICAL_H
    #define _DWS_NUMERICAL_H

    #include <limits>
    #include <iterator>
    
    namespace DWS
    {
    
    /*! 
     *\defgroup dws_numerical DWS Numerical Library
     * @{
     */
    
            // type definitions
        typedef double REAL_TYPE;  // can change to float
    
        // constants:
        const REAL_TYPE EPS(std::numeric_limits<REAL_TYPE>::epsilon());
        const REAL_TYPE REAL_MIN(std::numeric_limits<REAL_TYPE>::min());
        const REAL_TYPE PI(3.14159265358979323846);
        const REAL_TYPE E(2.718281828459045235360);
    
        // Library functions
        void SetMAXIT(int it);  //!< Sets maximum iterations for iterative functions.
    
        ////////////////////////////////////////////////////////////////////////
        // Special Mathematical functions
        ////////////////////////////////////////////////////////////////////////
        REAL_TYPE betai(REAL_TYPE a, REAL_TYPE b, REAL_TYPE x);//!< Incomplete Beta function (Numerical recipes p. 227)
        REAL_TYPE gammln(REAL_TYPE x); //!< Gamma fctn (Num. Rec. p. 214)
        REAL_TYPE gammp(REAL_TYPE a, REAL_TYPE x); //!< Incomplete Gamma Function (Num. Rec. p. 218)
        REAL_TYPE erff(REAL_TYPE x); //!< Error function (erf(x)) Num Rec. p 220
        REAL_TYPE factrl(int n); //!< Factorial: returns n! as floating pt num
        REAL_TYPE bico(int n, int k); //!< Binomial coefficient: returns binomial coefficient as floating pt num
        REAL_TYPE factln(int n); //!< Natural log: returns ln(n!)
        
    
        
        /*! \brief Probability density function evaluated at \a distance.
         *
         * Returns Gaussian probability function according to distance from mode
         * (\a distance: d) and standard deviation (\a sd: s).  Function returns
         * \f[ \frac{e^{\left( -d^2/2s^2\right)}}{\sqrt{2\pi}} \f]
         */
        REAL_TYPE normal_probability(REAL_TYPE distance, REAL_TYPE sd);
        
        /*! brief Round to significant digits.
         *
         *   Return x rounded to the specified number of significant digits, n, as
         *	counted from the first non-zero digit. 
         *	
         *	If n=0 (the default value for round2) then the magnitude of the
         *	number will be returned (e.g. round2(12) returns 10.0).  
         *	
         *	If n<0 then x will be rounded to the nearest multiple of n which, by 
         *	default, will be rounded to 1 digit (e.g. round2(1.23,-.28) will round 
         *	1.23 to the nearest multiple of 0.3.
         *
         *	Regardless of n, if x=0, 0 will be returned.
         */
        
///////////////////////////////////////////////////////////////////////////////
// Classes and structs
///////////////////////////////////////////////////////////////////////////////

        // Struct: QuadraticT
        // ---------------
        /*!
         * class for solving quadratic equations. the struct holds a, b and c 
         * parameters to a quadratic function.
         * The Solve() method fills \a first and \a second with the two 
         * solutions to the quadratic equation.
         */
        struct QuadraticT
        {
            REAL_TYPE _a, _b, _c;
            QuadraticT() {_a = _b = _c = 0;}
            
            /*! 
             * Fills \a first and \a second with two roots of quadratic
             * function
             */
            void Solve(REAL_TYPE& first, REAL_TYPE& second ) const;
            
            /*!
             * Returns derivative of the quadratic function.
             */
            REAL_TYPE DSolve() const;
        };


  /*@}*/
  
    }; /* namespace DWS */
    
    
#endif
// _DWS_MATH_H
