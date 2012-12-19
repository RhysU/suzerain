/***********************************************************************\
 *
 * File:           rngstream.cpp for multiple streams of Random Numbers
 * Language:       C++ (ISO 1998)
 * Copyright:      Pierre L'Ecuyer, University of Montreal
 * Notice:         This code can be used freely for personal, academic,
 *                 or non-commercial purposes. For commercial purposes,
 *                 please contact P. L'Ecuyer at: lecuyer@iro.umontreal.ca
 * Date:           14 August 2001
 *
 * $Id$
 *
\***********************************************************************/

/** @file
 * @copydoc rngstream.hpp
 */

#include <suzerain/rngstream.hpp>

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

namespace suzerain {

namespace
{
const double m1   =       4294967087.0;
const double m2   =       4294944443.0;
const double norm =       1.0 / (m1 + 1.0);
const double a12  =       1403580.0;
const double a13n =       810728.0;
const double a21  =       527612.0;
const double a23n =       1370589.0;
const double two17 =      131072.0;
const double two53 =      9007199254740992.0;
const double fact =       5.9604644775390625e-8;     /* 1 / 2^24  */

// The following are the transition matrices of the two MRG components
// (in matrix form), raised to the powers -1, 1, 2^76, and 2^127, resp.

const double InvA1[3][3] = {          // Inverse of A1p0
       { 184888585.0,   0.0,  1945170933.0 },
       {         1.0,   0.0,           0.0 },
       {         0.0,   1.0,           0.0 }
       };

const double InvA2[3][3] = {          // Inverse of A2p0
       {      0.0,  360363334.0,  4225571728.0 },
       {      1.0,          0.0,           0.0 },
       {      0.0,          1.0,           0.0 }
       };

const double A1p0[3][3] = {
       {       0.0,        1.0,       0.0 },
       {       0.0,        0.0,       1.0 },
       { -810728.0,  1403580.0,       0.0 }
       };

const double A2p0[3][3] = {
       {        0.0,        1.0,       0.0 },
       {        0.0,        0.0,       1.0 },
       { -1370589.0,        0.0,  527612.0 }
       };

const double A1p76[3][3] = {
       {      82758667.0, 1871391091.0, 4127413238.0 },
       {    3672831523.0,   69195019.0, 1871391091.0 },
       {    3672091415.0, 3528743235.0,   69195019.0 }
       };

const double A2p76[3][3] = {
       {    1511326704.0, 3759209742.0, 1610795712.0 },
       {    4292754251.0, 1511326704.0, 3889917532.0 },
       {    3859662829.0, 4292754251.0, 3708466080.0 }
       };

const double A1p127[3][3] = {
       {    2427906178.0, 3580155704.0,  949770784.0 },
       {     226153695.0, 1230515664.0, 3580155704.0 },
       {    1988835001.0,  986791581.0, 1230515664.0 }
       };

const double A2p127[3][3] = {
       {    1464411153.0,  277697599.0, 1610723613.0 },
       {      32183930.0, 1464411153.0, 1022607788.0 },
       {    2824425944.0,   32183930.0, 2093834863.0 }
       };



//-------------------------------------------------------------------------
// Return (a*s + c) MOD m; a, s, c and m must be < 2^35
//
double MultModM (double a, double s, double c, double m)
{
    double v;
    long a1;

    v = a * s + c;

    if (v >= two53 || v <= -two53) {
        a1 = static_cast<long> (a / two17);    a -= a1 * two17;
        v  = a1 * s;
        a1 = static_cast<long> (v / m);     v -= a1 * m;
        v = v * two17 + a * s + c;
    }

    a1 = static_cast<long> (v / m);
    /* in case v < 0)*/
    if ((v -= a1 * m) < 0.0) return v += m;   else return v;
}


//-------------------------------------------------------------------------
// Compute the vector v = A*s MOD m. Assume that -m < s[i] < m.
// Works also when v = s.
//
void MatVecModM (const double A[3][3], const double s[3], double v[3],
                 double m)
{
    int i;
    double x[3];               // Necessary if v = s

    for (i = 0; i < 3; ++i) {
        x[i] = MultModM (A[i][0], s[0], 0.0, m);
        x[i] = MultModM (A[i][1], s[1], x[i], m);
        x[i] = MultModM (A[i][2], s[2], x[i], m);
    }
    for (i = 0; i < 3; ++i)
        v[i] = x[i];
}


//-------------------------------------------------------------------------
// Compute the matrix C = A*B MOD m. Assume that -m < s[i] < m.
// Note: works also if A = C or B = C or A = B = C.
//
void MatMatModM (const double A[3][3], const double B[3][3],
                 double C[3][3], double m)
{
    int i, j;
    double V[3], W[3][3];

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j)
            V[j] = B[j][i];
        MatVecModM (A, V, V, m);
        for (j = 0; j < 3; ++j)
            W[j][i] = V[j];
    }
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            C[i][j] = W[i][j];
}


//-------------------------------------------------------------------------
// Compute the matrix B = (A^(2^e) Mod m);  works also if A = B.
//
void MatTwoPowModM (const double A[3][3], double B[3][3], double m, long e)
{
   int i, j;

   /* initialize: B = A */
   if (A != B) {
      for (i = 0; i < 3; ++i)
         for (j = 0; j < 3; ++j)
            B[i][j] = A[i][j];
   }
   /* Compute B = A^(2^e) mod m */
   for (i = 0; i < e; i++)
      MatMatModM (B, B, B, m);
}


//-------------------------------------------------------------------------
// Compute the matrix B = (A^n Mod m);  works even if A = B.
//
void MatPowModM (const double A[3][3], double B[3][3], double m, long n)
{
    int i, j;
    double W[3][3];

    /* initialize: W = A; B = I */
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j) {
            W[i][j] = A[i][j];
            B[i][j] = 0.0;
        }
    for (j = 0; j < 3; ++j)
        B[j][j] = 1.0;

    /* Compute B = A^n mod m using the binary decomposition of n */
    while (n > 0) {
        if (n % 2) MatMatModM (W, B, B, m);
        MatMatModM (W, W, W, m);
        n /= 2;
    }
}


//-------------------------------------------------------------------------
// Check that the seeds are legitimate values. Returns 0 if legal seeds,
// -1 otherwise.
//
int CheckSeed (const unsigned long seed[6])
{
    int i;

    for (i = 0; i < 3; ++i) {
        if (seed[i] >= m1) {
            cerr << "****************************************\n"
                 << "ERROR: Seed[" << i << "] >= 4294967087, Seed is not set."
                 << "\n****************************************\n\n";
            return (-1);
        }
    }
    for (i = 3; i < 6; ++i) {
        if (seed[i] >= m2) {
            cerr << "*****************************************\n"
                 << "ERROR: Seed[" << i << "] >= 4294944443, Seed is not set."
                 << "\n*****************************************\n\n";
            return (-1);
        }
    }
    if (seed[0] == 0 && seed[1] == 0 && seed[2] == 0) {
         cerr << "****************************\n"
              << "ERROR: First 3 seeds = 0.\n"
              << "****************************\n\n";
         return (-1);
    }
    if (seed[3] == 0 && seed[4] == 0 && seed[5] == 0) {
         cerr << "****************************\n"
              << "ERROR: Last 3 seeds = 0.\n"
              << "****************************\n\n";
         return (-1);
    }

    return 0;
}

} // end of anonymous namespace


//-------------------------------------------------------------------------
// Generate the next random number.
//
double rngstream::U01 ()
{
    long k;
    double p1, p2, u;

    /* Component 1 */
    p1 = a12 * Cg[1] - a13n * Cg[0];
    k = static_cast<long> (p1 / m1);
    p1 -= k * m1;
    if (p1 < 0.0) p1 += m1;
    Cg[0] = Cg[1]; Cg[1] = Cg[2]; Cg[2] = p1;

    /* Component 2 */
    p2 = a21 * Cg[5] - a23n * Cg[3];
    k = static_cast<long> (p2 / m2);
    p2 -= k * m2;
    if (p2 < 0.0) p2 += m2;
    Cg[3] = Cg[4]; Cg[4] = Cg[5]; Cg[5] = p2;

    /* Combination */
    u = ((p1 > p2) ? (p1 - p2) * norm : (p1 - p2 + m1) * norm);

    return (anti == false) ? u : (1 - u);
}


//-------------------------------------------------------------------------
// Generate the next random number with extended (53 bits) precision.
//
double rngstream::U01d ()
{
    double u;
    u = U01();
    if (anti) {
        // Don't forget that U01() returns 1 - u in the antithetic case
        u += (U01() - 1.0) * fact;
        return (u < 0.0) ? u + 1.0 : u;
    } else {
        u += U01() * fact;
        return (u < 1.0) ? u : (u - 1.0);
    }
}


//*************************************************************************
// Public members of the class start here


//-------------------------------------------------------------------------
// The default seed of the package; will be the seed of the first
// declared rngstream, unless SetPackageSeed is called.
//
double rngstream::nextSeed[6] =
{
   12345.0, 12345.0, 12345.0, 12345.0, 12345.0, 12345.0
};


//-------------------------------------------------------------------------
// constructor
//
rngstream::rngstream (const char *s) : name (s)
{
   anti = false;
   incPrec = false;

   /* Information on a stream. The arrays {Cg, Bg, Ig} contain the current
   state of the stream, the starting state of the current SubStream, and the
   starting state of the stream. This stream generates antithetic variates
   if anti = true. It also generates numbers with extended precision (53
   bits if machine follows IEEE 754 standard) if incPrec = true. nextSeed
   will be the seed of the next declared rngstream. */

   for (int i = 0; i < 6; ++i) {
      Bg[i] = Cg[i] = Ig[i] = nextSeed[i];
   }

   MatVecModM (A1p127, nextSeed, nextSeed, m1);
   MatVecModM (A2p127, &nextSeed[3], &nextSeed[3], m2);
}


//-------------------------------------------------------------------------
// Reset Stream to beginning of Stream.
//
void rngstream::ResetStartStream ()
{
   for (int i = 0; i < 6; ++i)
      Cg[i] = Bg[i] = Ig[i];
}


//-------------------------------------------------------------------------
// Reset Stream to beginning of SubStream.
//
void rngstream::ResetStartSubstream ()
{
   for (int i = 0; i < 6; ++i)
      Cg[i] = Bg[i];
}


//-------------------------------------------------------------------------
// Reset Stream to NextSubStream.
//
void rngstream::ResetNextSubstream ()
{
   MatVecModM(A1p76, Bg, Bg, m1);
   MatVecModM(A2p76, &Bg[3], &Bg[3], m2);
   for (int i = 0; i < 6; ++i)
       Cg[i] = Bg[i];
}


//-------------------------------------------------------------------------
bool rngstream::SetPackageSeed (const unsigned long seed[6])
{
   if (CheckSeed (seed))
      return false;                   // FAILURE
   for (int i = 0; i < 6; ++i)
#pragma warning(push,disable:2259)
      nextSeed[i] = seed[i];
#pragma warning(pop)
   return true;                       // SUCCESS
}


//-------------------------------------------------------------------------
bool rngstream::SetSeed (const unsigned long seed[6])
{
   if (CheckSeed (seed))
      return false;                   // FAILURE
   for (int i = 0; i < 6; ++i)
#pragma warning(push,disable:2259)
      Cg[i] = Bg[i] = Ig[i] = seed[i];
#pragma warning(pop)
   return true;                       // SUCCESS
}


//-------------------------------------------------------------------------
// if e > 0, let n = 2^e + c;
// if e < 0, let n = -2^(-e) + c;
// if e = 0, let n = c.
// Jump n steps forward if n > 0, backwards if n < 0.
//
void rngstream::AdvanceState (long e, long c)
{
    double B1[3][3], C1[3][3], B2[3][3], C2[3][3];

    if (e > 0) {
        MatTwoPowModM (A1p0, B1, m1, e);
        MatTwoPowModM (A2p0, B2, m2, e);
    } else if (e < 0) {
        MatTwoPowModM (InvA1, B1, m1, -e);
        MatTwoPowModM (InvA2, B2, m2, -e);
    }

    if (c >= 0) {
        MatPowModM (A1p0, C1, m1, c);
        MatPowModM (A2p0, C2, m2, c);
    } else {
        MatPowModM (InvA1, C1, m1, -c);
        MatPowModM (InvA2, C2, m2, -c);
    }

    if (e) {
        MatMatModM (B1, C1, C1, m1);
        MatMatModM (B2, C2, C2, m2);
    }

    MatVecModM (C1, Cg, Cg, m1);
    MatVecModM (C2, &Cg[3], &Cg[3], m2);
}


//-------------------------------------------------------------------------
void rngstream::GetState (unsigned long seed[6]) const
{
   for (int i = 0; i < 6; ++i)
      seed[i] = static_cast<unsigned long> (Cg[i]);
}


//-------------------------------------------------------------------------
void rngstream::WriteState (std::ostream &out) const
{
    out << "The current state of the Rngstream";
    if (name.size() > 0)
        out << " " << name;
    out << ":\n   Cg = { ";

    for (int i = 0; i < 5; i++) {
        out << static_cast<unsigned long> (Cg [i]) << ", ";
    }
    out << static_cast<unsigned long> (Cg [5]) << " }\n\n";
}


//-------------------------------------------------------------------------
void rngstream::WriteStateFull (std::ostream &out) const
{
    int i;

    out << "The RngStream";
    if (name.size() > 0)
        out << " " << name;
    out << ":\n   anti = " << (anti ? "true" : "false") << "\n";
    out << "   incPrec = " << (incPrec ? "true" : "false") << "\n";

    out << "   Ig = { ";
    for (i = 0; i < 5; i++) {
        out << static_cast<unsigned long> (Ig [i]) << ", ";
    }
    out << static_cast<unsigned long> (Ig [5]) << " }\n";

    out << "   Bg = { ";
    for (i = 0; i < 5; i++) {
        out << static_cast<unsigned long> (Bg [i]) << ", ";
    }
    out << static_cast<unsigned long> (Bg [5]) << " }\n";

    out << "   Cg = { ";
    for (i = 0; i < 5; i++) {
        out << static_cast<unsigned long> (Cg [i]) << ", ";
    }
    out << static_cast<unsigned long> (Cg [5]) << " }\n\n";
}


//-------------------------------------------------------------------------
void rngstream::IncreasedPrecis (bool incp)
{
   incPrec = incp;
}


//-------------------------------------------------------------------------
void rngstream::SetAntithetic (bool a)
{
   anti = a;
}


//-------------------------------------------------------------------------
// Generate the next random number.
//
double rngstream::RandU01 ()
{
   if (incPrec)
      return U01d();
   else
      return U01();
}


//-------------------------------------------------------------------------
// Generate the next random integer.
//
int rngstream::RandInt (int low, int high)
{
#pragma warning(push,disable:2259)
    return low + static_cast<int> ((high - low + 1.0) * RandU01 ());
#pragma warning(pop)
}

//-------------------------------------------------------------------------
// Approximate the inverse of the standard normal per Peter John Acklam
// (http://home.online.no/~pjacklam/notes/invnorm/) with error handling
// following Chad Sprouse's implementation given at
// (http://home.online.no/~pjacklam/notes/invnorm/impl/sprouse/ltqnorm.c).
//
static double acklam(const double p)
{
   // Coefficients in rational approximations.
   static const double a1 = -3.969683028665376e+01;
   static const double a2 =  2.209460984245205e+02;
   static const double a3 = -2.759285104469687e+02;
   static const double a4 =  1.383577518672690e+02;
   static const double a5 = -3.066479806614716e+01;
   static const double a6 =  2.506628277459239e+00;

   static const double b1 = -5.447609879822406e+01;
   static const double b2 =  1.615858368580409e+02;
   static const double b3 = -1.556989798598866e+02;
   static const double b4 =  6.680131188771972e+01;
   static const double b5 = -1.328068155288572e+01;

   static const double c1 = -7.784894002430293e-03;
   static const double c2 = -3.223964580411365e-01;
   static const double c3 = -2.400758277161838e+00;
   static const double c4 = -2.549732539343734e+00;
   static const double c5 =  4.374664141464968e+00;
   static const double c6 =  2.938163982698783e+00;

   static const double d1 =  7.784695709041462e-03;
   static const double d2 =  3.224671290700398e-01;
   static const double d3 =  2.445134137142996e+00;
   static const double d4 =  3.754408661907416e+00;

   // Define break-points.
   const double p_low  = 0.02425;
   const double p_high = 1 - p_low;

   // Compute rational approximation for...
   double x, q, r;
   if (p < 0) {                      // ...domain error.
      errno = EDOM;
      return 0;
#pragma warning(push,disable:1572)
   } else if (p == 0) {              // ...improbable point.
#pragma warning(pop)
      errno = ERANGE;
#ifdef HUGE_VAL
      return -HUGE_VAL;
#else
      return -DBL_MAX;
#endif
   } else if (p < p_low) {           // ...lower region.
      q = sqrt(-2*log(p));
      x = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) /
           ((((d1*q + d2)*q + d3)*q + d4)*q + 1);
   } else if (p <= p_high) {         // ...central region
      q = p - 0.5;
      r = q*q;
      x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6)*q /
          (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1);
   } else if (p < 1) {               // ...upper region.
      q = sqrt(-2*log(1 - p));
      x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) /
            ((((d1*q + d2)*q + d3)*q + d4)*q + 1);
#pragma warning(push,disable:1572)
   } else if (p == 1) {              // ...improbable point.
#pragma warning(pop)
      errno = ERANGE;
#ifdef HUGE_VAL
      return HUGE_VAL;
#else
      return DBL_MAX;
#endif
   } else {                          // ...domain error.
      errno = EDOM;
      return 0;
   }

   // One iteration of Halley's rational method improves to machine precision.
   static const double neg_1_sqrt2 = - 1 / M_SQRT2;
   static const double sqrt_2pi    = sqrt(2 * M_PI);
   double e, u;
   e = 0.5 * erfc(x * neg_1_sqrt2) - p;
   u = e * sqrt_2pi * exp(x*x/2);
   x = x - u/(1 + x*u/2);

   return x;
}

//-------------------------------------------------------------------------
// Generate the next random number from the standard normal
//
double rngstream::RandN01 ()
{
   return acklam(RandU01());
}

} // end namespace suzerain
