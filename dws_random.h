/******************************************************************************/
// FILE NAME: dws_random.h
// VERSION: $Revision: 1.2 $
// REMARKS: Random number library
//------------------------------------------------------------------------------
// Copyright © 2002 - 2005 
//------------------------------------------------------------------------------
// AUTHOR/S:  Dylan W. Schwilk 
//------------------------------------------------------------------------------
// CREATION DATE: 020830
//---------+------+---------------------------------------------+---------------
// REVISION LOG:
// $Log: dws_random.h $
// Revision 1.2  2005/01/26 17:34:50  schwilkd
// Added Doxygen comments. Renamed functions to use lower_case naming
// style to match all my other libraries.
//
// Revision 1.1  2005/01/26 17:34:14  schwilkd
// Initial revision
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
 *
 * Library with random number functions. The functions RandomReal and
 * Random integer return random deviates from a uniform distribution
 * in the range requested.  Use the Randomize() function to seed the
 * number generator from the system clock.
 *
 * This library treats random numbers from a functional, not object,
 * perspective.  Static data members are saved at the module level and
 * this library may not be the best option if code needs alternating
 * access to many different types of generators.  In such a case, a
 * random number class library like that available at www.boost.org
 * may may more sense.
 */
 
#ifndef DWS_RANDOM_H
    #define DWS_RANDOM_H
    
    #include "dws_numerical.h"

    namespace DWS
    {   
    
    /*! \defgroup dws_random DWS Random Number Library
     * @{
     */
 
        /* Function: raw_rand()
         * -------------------*/
         /*!
         * Returns a random number between 0 and 1 exclusive
         */
        inline REAL_TYPE raw_rand();
    
        /* Function: exp_deviate();
         * -----------------------*/
        /*! This function returns a exponentially distributed, random
         * deviate (x) with mean 1/lambda, The probability density
         * function from which the result is drawn follows the
         * formula: \f$ P(x) = \frac{1}{\lambda} e^\frac{-\lambda}{x}
         * \f$ .
         *
         */
        REAL_TYPE exp_deviate(REAL_TYPE lambda);
    
        /* Function normal_deviate(mean, variance)
         * --------------------------------------*/
        /*!
         * This function returns a normal deviate from a distribution
         * with the specified mean and variance
         */
        REAL_TYPE normal_deviate(REAL_TYPE mean, REAL_TYPE variance);
        
        /* Function poission_deviate(REAL_TYPE xm)
         * --------------------------------------*/
        /*! 
         * Returns integer value (as floating point #) of random poisson
         * deviate with mean \a xm.
         */
        REAL_TYPE poission_deviate(REAL_TYPE xm);
    
        /* Function: randomize();
         * --------------------------------------*/
        /*! 
         * Sets the random number seed from clock. Returns seed.
         */
        long randomize(void);
    
        /* Function: set_randomizer_seed(int s);
         * --------------------------------------*/
        /*!
         * Sets the random number seed
         */
        void set_randomizer_seed(int s);
    
    
        /* Function: random_angle_sin_cos(outSin, outCos)
         * --------------------------------------*/
        /*!
         * This function returns the sine and cosine of a random angle.  It is
         * more efficient than calling sin() and cos() on a uniform random number
         * between 0 and 2pi.
         */
        void random_angle_sin_cos(REAL_TYPE & outSin, REAL_TYPE & outCos);
    
    
        /* Function: random_bit()
         * --------------------------------------*/
       /*! Returns an integer 0 or 1 bernouli .5
         */
        int random_bit();
    
        /* Function: random_integer
         * Usage: n = random_integer(low, high);
         * --------------------------------------*/
        /*! 
         * Returns a random integer in the range low to high, inclusive.
         */
        int random_integer(int low, int high);
    
        /* Function: random_real(low, high);
         * --------------------------------------*/
        /*! 
         *Returns a random real number between low and high (exclusive of the
         * endpoint values).
         */
        double random_real(REAL_TYPE low, REAL_TYPE high);
    

    
        /* Function: random_chance
         * Usage if (random_chance(p)) . . .
         * --------------------------------------*/
        /*!
         * returns TRUE with probability indicated by p
         */
        bool random_chance(REAL_TYPE p);
    
    
        // Class TIRand
        // ------------
        /*! Integer random number generator fxn object.
         *
         * Function object for use in STL algorithms that require a random
         * number generator.  Generates random integers in range 0 - \a n.
         */
        class TIRand {
         public :
         inline int operator()(int n) {return random_integer(0,n-1);};
        };   
        
    /*@}*/
    
    
    }; // end namespace DWS
#endif  // DWS_RANDOM_H
