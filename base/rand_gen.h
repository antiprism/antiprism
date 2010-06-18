/* Public Domain
 * Robin Whittle  rw@firstpr.com.au    2005 September 21 
 * http://www.firstpr.com.au/dsp/rand31/rand31-park-miller-carta.cc.txt
 *
 * Added to Antiprism - http://www.antiprism.com by Adrian Rossiter.
 * The original commented code is included at the end of this file.
 */

/// @cond ALL

/*!\file rand_gen.h
 * \brief A pseudo-random number generator
 */

#ifndef RAND_GEN_H
#define RAND_GEN_H

#include <stdlib.h>
#include <math.h>
#include <time.h>

#define consta 16807            

///Psuedo-random number generator
class rand_gen {
   private:
      long unsigned int seed31;   
      long unsigned int nextrand();

   public:
      ///Constructor
      rand_gen() {seed31 = 1;}

      ///Constructor
      /**\param seedin the seed, if zero is used, then
       * the seed will be set to 1. */
      rand_gen(long unsigned int seedin) { seedi(seedin); }

      ///Set the seed.
      /**\param seedin the seed, if zero is used, then
       * the seed will be set to 1. */
      void seedi(long unsigned int seedin)
         { if (seedin == 0) seedin = 1; seed31 = seedin; }

      ///Set the seed with the current time.
      /**\param seedin the seed, if zero is used, then
       * the seed will be set to 1. */
      void time_seed()
         { seedi(time(0)); nextrand(); }

      ///Get a psuedo-random integer
      /**\return a pseudo-random integer. */
      long unsigned int ranlui() { return nextrand(); }

      ///Get a psuedo-random floating point number in range 0.0 <= num <= 1.0
      /**\return a pseudo-random float in the range 0.0 - 1.0, including 1.0.*/
      double ranf() { return (nextrand() / 2147483647.0); }
    
      ///Get a psuedo-random floating point number in range 0.0 <= num < 1.0
      /**\return a pseudo-random float in the range 0.0 - 1.0, excluding 1.0.*/
      double ranf_exclude_end() { return (nextrand() / 2147483648.0); }
    
      ///Get a psuedo-random floating point number in range low <= num <= high
      /**\return a pseudo-random float in the range low - high, including high.*/
      double ran_in_range(double low, double high)
      {
         return low + (high-low)*ranf();
      }
      
      ///Get a psuedo-random floating point number in range low <= num < high
      /**\return a pseudo-random float in the range low - high, excluding high.*/
      double ran_in_range_exclude_end(double low, double high)
      {
         return low + (high-low)*ranf_exclude_end();
      }

      ///Get a psuedo-random integer in range low <= num <=high
      /**\return a pseudo-random integer in the range low - high, including high.*/
      long ran_int_in_range(long low, long high)
      {
         return (long)floor(low + ((high+1)-low)*ranf_exclude_end());
      }
      
      ///Get a psuedo-random integer in range 0 <= num <high
      /**\return a pseudo-random integer in the range low - high, including high.*/
      long operator() (long high)
      {
         return (long)floor(high*ranf_exclude_end());
      }


};


inline long unsigned int rand_gen::nextrand()
{
   long unsigned int hi, lo;
   lo = consta * (seed31 & 0xFFFF);
   hi = consta * (seed31 >> 16);
   lo += (hi & 0x7FFF) << 16;
   lo += hi >> 15;
   if (lo > 0x7FFFFFFF) lo -= 0x7FFFFFFF;          
   return ( seed31 = (long)lo );       
}


/*
class rand31dc {
                                    // The sole item of state - a 32 bit 
                                    // integer.
    long unsigned int seed31;   

public:
                                    // Constructor sets seed31 to 1, in case 
                                    // no seedi operation is used.
    rand31dc() {seed31 = 1;}                                    
                                    
                                    // Declare methods.
                                    
    void              seedi(long unsigned int);
    long unsigned int ranlui(void);  
    float             ranf(void);
    

private:
                                    // nextrand()
                                    //
                                    // Generate next pseudo-random number.
                                    
                                    // Multiplier constant = 16807 = 7^5.  
                                    // This is 15 bits.
                                    //
                                    // Park and Miller in 1993 recommend
                                    // 48271, which they say produces a 
                                    // somewhat better quality of 
                                    // pseudo-random results.  
                                    //
                                    // 48271 can't be used with the 
                                    // following implementation of Carta's 
                                    // algorithm, since it is 16 bits and 
                                    // would result in bit 31 potentially 
                                    // being set in lo in the first
                                    // multiplication.  (A similar problem
                                    // occurs later with the higher bits of
                                    // hi.)

    #define consta 16807            
                                    // Modulus constant = 2^31 - 1 =
                                    // 0x7FFFFFFF.   We use this explicitly
                                    // in the code, rather than define it
                                    // somewhere, because this is a value
                                    // which must not be changed and should
                                    // always be recognised as a zero 
                                    // followed by 31 ones.
                                            
    long unsigned int nextrand()
    {
                                    // Two 32 bit integers for holding
                                    // parts of the (seed31 * consta)
                                    // multiplication product which would 
                                    // normally need a 46 bit word. 
                                    // 
                                    // lo 31 bits       30  -  0 
                                    // hi 30 bits   45  -  16  
                                    //
                                    // These overlap in their value.
                                    //
                                    // Bit 0 of hi has the same weight in 
                                    // the result as bit 16 of lo.
                                    
        long unsigned int hi, lo;

                                    // lo = 31 bits:
                                    //  
                                    //    low 16 bits (15-0) of seed31 
                                    //  * 15 bit consta 
                                     
        lo = consta * (seed31 & 0xFFFF);
        
                                    // hi = 30 bits:
                                    //
                                    //    high 15 bits (30-16) of seed31
                                    //  * 15 bit consta 
                                    
        hi = consta * (seed31 >> 16);
        
                                    // The new pseudo-random number is the 
                                    // 46 bit product mod(0x7FFFFFF).  Our
                                    // task is to calculate it with 32
                                    // bit words and maths, and without
                                    // division.
                                    //
                                    // The first section is easy to
                                    // understand.  We have a bunch of
                                    // bits - bits 14 to 0 of hi - 
                                    // which overlap with and carry the
                                    // same weight as bits 30 to 16 of
                                    // lo.
                                    //
                                    // Add the low 15 bits of hi into
                                    // bits 30-16 of lo.  
                                    
        lo += (hi & 0x7FFF) << 16;
        
                                    // The result may set bit 31 of lo, but
                                    // it will not overflow lo.
                                    //
                                    // So far, we got some of our total
                                    // result in lo.
                                    //
                                    // The only other part of the result
                                    // we need to deal with is bits
                                    // 29 to 15 of hi. 
                                    //
                                    // These bits carry weights of bits
                                    // 45 to 31 in the value of the 
                                    // multiplication product of the usual
                                    // Park-Miller algorithm.
                                    //
                                    // David Carta writes that in order
                                    // to get the mod(0x7FFFFFF) of the
                                    // multiplication product we should
                                    // simply add these bits into the
                                    // bit positions 14 to 0 of lo.
                                    //
        lo += hi >> 15;             //
                                    // In order to be able to get away with
                                    // this, and to perform the following
                                    // simple mod(0x7FFFFFFF) operation,
                                    // we need to be sure that the result 
                                    // of the addition will not exceed:
                                    //                                      
                                    // 2 * 0x7FFFFFFF = 0xFFFFFFFE
                                    // 
                                    // This is assured as per the diagrams
                                    // above.
                                    // Note that in the vast majority of 
                                    // cases, lo will be less than 
                                    // 0x7FFFFFFFF. 
                                    
        if (lo > 0x7FFFFFFF) lo -= 0x7FFFFFFF;          
        
                                    // lo contains the new pseudo-random
                                    // number.  Store it to the seed31 and
                                    // return it.
        
        return ( seed31 = (long)lo );       
    }
};

                                    /////////////////////////////////////////
                                    //
                                    // Implementations of the methods.
                                    
                                    // seedi()
                                    //
                                    // Set the seed from a long unsigned 
                                    // integer.  If zero is used, then
                                    // the seed will be set to 1.
                                                                        
void rand31dc::seedi(long unsigned int seedin)
{
    if (seedin == 0) seedin = 1;
    seed31 = seedin;
}
                                    
                                    // ranlui()
                                    //
                                    // Return next pseudo-random value as
                                    // a long unsigned integer.
                                    
long unsigned int rand31dc::ranlui(void)  
{
    return nextrand();
}

                                    // ranf()
                                    //
                                    // Return next pseudo-random value as
                                    // a floating point value.


float rand31dc::ranf(void)  
{
    return (nextrand() / 2147483647.0);
}


//////////////////////////////////////////////////////////////////////////////

*/


/// @endcond



#endif // RAND_GEN_H

