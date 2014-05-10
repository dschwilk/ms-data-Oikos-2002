/**************************************************************************/
// FILE NAME: dws_random.cpp
// VERSION: $Revision: 1.1 $
// REMARKS: Random number library. Implements dws_random.h   
//------------------------------------------------------------------------*/
// Copyright © 2000 - 2002 */
/*------------------------------------------------------------------------*/
// AUTHOR/S:  Dylan W. Schwilk */
//------------------------------------------------------------------------*/
// CREATION DATE: 20/09/02 16:54
//---------+---------+----------------------------------------------------*/
// REVISION LOG:
// $Log: dws_random.cpp $
// Revision 1.1  2005/01/26 17:39:12  schwilkd
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

#include "dws_random.h"

#include <ctime>
#include <cmath>
#include "dws_numerical.h"

namespace DWS
{
    //  Constants
  
    //  Static data members:
    static long seed = -1;   // random number seed set by Randomize()
    static long iseed = -1; // used by RandomBit()
    static bool GaussInReserve = false; // Flag that signals existence of reserve normal deviate
    static REAL_TYPE GaussSet = 0;   // reserve normal deviate used by GaussDev()


    /* Static functions
     */
    inline static REAL_TYPE Ran1(long *idum);
    static REAL_TYPE GaussDev(void);


    /* Implementation of public functions
     * ----------------
     */

    /* Function: SetRandomizerSeed(int s);
     * ----------------------
     */
    void set_randomizer_seed(int s)
    {
        iseed = seed = (s>0)? -s:s;
    }

    long randomize(void)
    {
        seed = ((long) std::time(NULL));
        if (seed > 0) seed = -seed;

        return iseed=seed;  // set iseed for RandomBit
    }

    /* RawRand
     * -------r
     */
    inline REAL_TYPE raw_rand()
    {
        return Ran1(&seed);
    }

    /* Function: RandomAngleSinCos(outSin, outCos)
     * ------------------------------
     */
    void random_angle_sin_cos(REAL_TYPE & outSin, REAL_TYPE & outCos)
    {
        REAL_TYPE rsq, v1, v2;

        do
        {      // get random deviates between -1 and 1
            v1 = REAL_TYPE(2.0*raw_rand()-1.0);
            v2 = REAL_TYPE(2.0*raw_rand()-1.0);
            rsq = v1*v1+v2*v2;  // check if they are in unit circle
        }
        while (rsq >= 1.0 || rsq == 0.0);

        outSin = REAL_TYPE(v1/std::sqrt(rsq));
        outCos = REAL_TYPE(v2/std::sqrt(rsq));
    }

    /* RandomBit()
     * -----------
     * Returns a random 0 or 1 each with probablility of 0.5.  Faster than using RandomChance.
     */
    int random_bit()
    {
    // constants for RandomBit function
    static const int IB1(1);
    static const int IB2(2);
    static const int IB5(16);
    static const int IB18(131072);
    static const int MASK(IB1+IB2+IB5);
    
        if (iseed & IB18)
        {// change all masked bits, shift and put 1 in bit 1
            iseed=((iseed ^ MASK) << 1) | IB1;
            return 1;
        }
        else
        {
            iseed <<= 1;
            return 0;
        }
    }


    /* RandomInteger
     * -------------
     */
    int random_integer(int low, int high)
    {
        int k;
        double d;

        d = raw_rand();
        k = (int) (d * ((high - low) + 1));
        return (low+k);
    }

    /* random_real
     * ----------
     */
    double random_real(REAL_TYPE low, REAL_TYPE high)
    {
        double d;

        d = raw_rand();
        d *= (high - low);
        return (low + d);
    }


    /* random_chance
     * ------------
     */
    bool random_chance(REAL_TYPE p)
    {
        return (raw_rand() < p);
    }

    /* Function: exp+deviate();
     * -----------------------
     * This function returns a n exponentially distributed, random deviate with mean 1/lambda,
     * using Ran1() as the source of the uniform distribution.
     */
    REAL_TYPE exp_deviate(REAL_TYPE lambda)
    {
        REAL_TYPE dum;

        do
        {
            dum=raw_rand();
        }
        while (dum == 0.0);
        return ( - std::log(dum) / lambda);
    }


    /* Function normal_deviate(mean, variance)
     * --------------------------------------
     */

    REAL_TYPE normal_deviate(REAL_TYPE mean, REAL_TYPE variance)
    {
        return (GaussDev()*variance + mean);
    }


    /* Function poission_deviate(REAL_TYPE xm)
     * --------------------------------------
     */
REAL_TYPE poission_deviate(REAL_TYPE xm)
{    
     static REAL_TYPE sq, alxm, g;
     static REAL_TYPE oldm = -1.0;  // flag to indicate if xm is different 
                                    // than last call
                                    // this is an obvious candidate for a
                                    // class 
    REAL_TYPE em,t,y;
                                    
     if(xm<12.0) {
        if(xm!=oldm) {
            oldm=xm;
            g=std::exp(-xm);
        }
        em=-1;
        t=1.0;
        do {
            ++em;
            t *= raw_rand();
        } while(t>g);
    } else {
        if(xm!=oldm) {
            oldm=xm;
            sq=std::sqrt(2.0*xm);
            alxm=log(xm);
            g=xm*alxm-gammln(xm+1.0);
        }
        do {
            do {
                y=tan(DWS::PI*raw_rand());
                em=sq*y+xm;
            } while(em<0.0);
            em=std::floor(em);
            t=0.9*(1.0+y*y)*std::exp(em*alxm-gammln(em+1.0)-g);
        }while(raw_rand() > t);
    }
    return em;
}

    /*----------------------------------------------------------------*/
    /* Private (static) functions */


    /* Function: Ran1
     * Usage: x= Ran1(*idum);
     * ---------------------
     * "Minimal" random number generator of Park and Miller with Bays=Durham
     * shuffle and added safeguards.  Returns a uniform random deviate between
     * 0 and 1.0 exclusive.  Call with idum a negative integer to initialize,
     * thereafter to not alter idum between successive deviates in a sequence.
     */
    REAL_TYPE Ran1(long *idum)
    {
    const long IA(16807);
    const long IM(2147483647);
    const REAL_TYPE  RNMX(1.0-EPS);
    const REAL_TYPE AM(1.0/IM);
    const long  IQ(127773);
    const long IR(2836);
    const long NTAB(32);
    const long  NDIV(1+(IM-1)/NTAB);

    
        int j;
        long k;
        static long iy = 0;
        static long iv[NTAB];
        REAL_TYPE temp;

        if (*idum <=0 || !iy)
        {  // initialize
            if (-(*idum) < 1) *idum = 1;
            else *idum = -(*idum);
            for (j=NTAB+7; j>=0; j--)
            {
                k=(*idum)/IQ;
                *idum=IA*(*idum-k*IQ)-IR*k;
                if (*idum<0) *idum += IM;
                if (j < NTAB) iv[j] = *idum;
            }
            iy=iv[0];
        }
        k=(*idum)/IQ;     // Start here when not initializing
        *idum=IA*(*idum-k*IQ)-IR*k;  // Compute idum=(IA*idum) % IM without
        if (*idum < 0) *idum +=IM;  // overflows by Schrage's method.
        j=iy/NDIV;      // Will be in the range 0..NTAB-1
        iy=iv[j];
        iv[j] = *idum;
        if ((temp = (REAL_TYPE) AM*iy) > RNMX) return (REAL_TYPE)(RNMX);
        else return (temp);
    }


    /* Function: GaussDev
     * -------------------
     * This function returns a normally distributed deviate w/ zero mean and
     * unit variance, using Ran1() as the source of the uniform deviates.
     * Since the function creates two normal deviates each call, the extra is
     * held in reserve in static variable GaussSet (signalled by bool
     * GaussInReserve == true) to improve efficiency.
     */

    REAL_TYPE GaussDev(void)
    {
        REAL_TYPE fac, rsq, v1, v2;

        if (!GaussInReserve)
        {
            do
            {
                v1 = REAL_TYPE(2.0*raw_rand()-1.0);
                v2 = REAL_TYPE(2.0*raw_rand()-1.0);
                rsq = v1*v1+v2*v2;  // check if they are in unit circle
            }
            while (rsq >= 1.0 || rsq == 0.0);
            fac = REAL_TYPE(std::sqrt(-2.0*std::log(rsq)/rsq));
            // Box-Muller transformation to get two normal deviates.  Return
            // one and save other for next function call.
            GaussSet = v1*fac;
            GaussInReserve = true;  // Set flag
            return v2*fac;
        }
        else
        {
            GaussInReserve = false;  // unset flag
            return GaussSet;
        }
    }

}; // namespace DWS
