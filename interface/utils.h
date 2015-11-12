#ifndef _INC_UTILS
#define _INC_UTILS

#include <vector>
#define _USE_MATH_DEFINES // include M_PI constant
#include <cmath>

#define SQ(x) (x*x)
#define CB(x) (x*x*x)
#define QU(x) (x*x*x*x)

//#define M_PI 3.1415927410125732421875

// Set option flags for CUBA integrator
// 
// smoothing only used by Suave
// smoothing and takeOnlyGridFromFile only used by Vegas
//
// verbosity = 0-3
// subregion true = only last set of samples is used for final evaluation of integral
// smoothing true = smoothe importance function
// retainStateFile true => retain state file when integration ends
// takeOnlyGridFromFile false => full state taken from file (if present), true => only grid is taken (e.g. to use it for another integrand)
// level determines random-number generator:
//  seed = 0 => Sobol quasi-random
//  seed != 0 => level is used:
//    level = 0 => Mersenne Twister pseudo-random
//    level != 0 => Ranlux pseudo-random, level = determines period of generator
unsigned int setFlags(char verbosity = 0, bool subregion = false, bool retainStateFile = false, unsigned int level = 0, bool smoothing = false, bool takeOnlyGridFromFile = true);

template<typename T> T sign(const T x){
  if(x > 0)
    return 1;
  else if(!x)
    return 0;
  else
    return -1;
}

// Used to compute Jacobian for Transfer Function
inline double dEoverdP(const double E, const double m){
  const double rad = SQ(E) - SQ(m);
  if(rad <= 0)
    return 0.;
  else
    return E/sqrt(rad);
}

// We flatten the Breit-Wigners by doing a change of variable for each resonance separately
// s = M G tan(y) + M^2
// jac = M G / cos^2(y)
// ==> BW(s(y))*jac(y) is flat in the variable y, as BW(s) = 1/((s-M^2)^2 - (GM)^2)
// Where y = -arctan(M/G) + (pi/2+arctan(M/G))*x_foam (x_foam between 0 and 1 => s between 0 and infinity)
inline void flattenBW(const double psPoint, const double mass, const double width, double &s, double &jac){
  const double range = M_PI/2. + atan(mass / width);
  const double y = - atan(mass / width) + range * psPoint;
  s = mass * width * tan(y) + SQ(mass);
  jac = range * mass * width / SQ(cos(y));
}

// If the NWA is used for a particle, a multiplication factor has to be introduced
// because of the integrated-out delta function
inline double jacobianNWA(const double mass, const double width){
  return ( M_PI/2. + atan(mass / width) ) * mass * width;
}

// Compute cos(x +- 2*pi/3) in a more "analytical" way (pm = +- 1)
// Useful for solveCubic
inline double cosXpm2PI3(const double x, const double pm){
  return -0.5*( cos(x) + pm * sin(x) * sqrt(3.) );
}

// Finds the real solutions to a*x^2 + b*x + c = 0
// Uses a numerically more stable way than the "classroom" method.
// Handles special cases a=0 and/or b=0.
// Appends the solutions to the std::vector roots, making no attempt to check whether the vector is empty.
// Double roots are present twice.
// If verbose is true (default is false), the solutions are printed, 
// as well as the polynomial evaluated on these solutions
//
// See https://fr.wikipedia.org/wiki/Équation_du_second_degré#Calcul_numérique
bool solveQuadratic(const double a, const double b, const double c, 
                    std::vector<double>& roots, 
                    bool verbose = false
                    );

// Finds the real solutions to a*x^3 + b*x^2 + c*x + d = 0
// Handles special case a=0.
// Appends the solutions to the std::vector roots, making no attempt to check whether the vector is empty.
// Multiple roots appear multiple times.
// If verbose is true (default is false), the solutions are printed, 
// as well as the polynomial evaluated on these solutions
//
// Inspired by "Numerical Recipes" (Press, Teukolsky, Vetterling, Flannery), 2007 Cambridge University Press
bool solveCubic(const double a, const double b, const double c, const double d, 
                std::vector<double>& roots, 
                bool verbose = false
                );

// Finds the real solutions to a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
// Handles special case a=0.
// Appends the solutions to the std::vector roots, making no attempt to check whether the vector is empty.
// Multiple roots appear multiple times.
// If verbose is true (default is false), the solutions are printed, 
// as well as the polynomial evaluated on these solutions
//
// See https://en.wikipedia.org/wiki/Quartic_function#Solving_by_factoring_into_quadratics
//    https://fr.wikipedia.org/wiki/Méthode_de_Descartes
//
// The idea is to make a change of variable to eliminate the term in x^3, which gives a "depressed quartic",
// then to try and factorize this quartic into two quadratic equations, which each then gives up to two roots.
// The factorization relies on solving a cubic equation (which is always possible), 
// then taking the square root of one of these solution (ie there must be a positive solution).
bool solveQuartic(const double a, const double b, const double c, const double d, const double e, 
                  std::vector<double>& roots, 
                  bool verbose = false
                  );

// Solves the system:
// a20*E1^2 + a02*E2^2 + a11*E1*E2 + a10*E1 + a01*E2 + a00 = 0
// b20*E1^2 + b02*E2^2 + b11*E1*E2 + b10*E1 + b01*E2 + b00 = 0
// Which corresponds to finding the intersection points of two conics.
// Appends the (x,y) solutions to the std::vectors E1, E2, making no attempt to check 
// whether these vectors are empty.
// In most cases it simply comes down to solving a quartic equation:
//   - eliminate the E1^2 term
//   - solve for E1 (linear!)
//   - inserting back gives a quartic function of E2
//   - find solutions for E1
// The procedure becomes tricky in some special cases
// (intersections aligned along x- or y-axis, degenerate conics, ...)
bool solve2Quads(const double a20, const double a02, const double a11, const double a10, const double a01, const double a00,
                 const double b20, const double b02, const double b11, const double b10, const double b01, const double b00,
                 std::vector<double>& E1, std::vector<double>& E2, 
                 bool verbose = false
                 );

// Solves the system:
// a11*E1*E2 + a10*E1 + a01*E2 + a00 = 0
// b11*E1*E2 + b10*E1 + b01*E2 + b00 = 0
// Which corresponds to finding the intersection points of two conics.
// Appends the (x,y) solutions to the std::vectors E1, E2, making no attempt to check 
// whether these vectors are empty.
bool solve2QuadsDeg(const double a11, const double a10, const double a01, const double a00,
                    const double b11, const double b10, const double b01, const double b00,
                    std::vector<double>& E1, std::vector<double>& E2, 
                    bool verbose = false
                    );

// Solves the system:
// a10*E1 + a01*E2 + a00 = 0
// b10*E1 + b01*E2 + b00 = 0
// Appends the (x,y) solutions to the std::vectors E1, E2, making no attempt to check 
// whether these vectors are empty.
bool solve2Linear(const double a10, const double a01, const double a00,
                  const double b10, const double b01, const double b00,
                  std::vector<double>& E1, std::vector<double>& E2, 
                  bool verbose = false);


double BreitWigner(const double s, const double m, const double g);

#endif
