#ifndef SPH_GLOBAL_H_
#define SPH_GLOBAL_H_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;

#define real			double
#define OK				0
#define ERROR			1
#define PI 				3.141592654
#define LINKEDLISTEND	-1000000
#define ROOT_PROCESS	0
#define NOPROC			-100000
#define LENNARDJONES	-1

// Useful macros
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define AVG(x,y) ( 0.5 * ( (x)+(y) ) )
#define SGN(x) ( (x>=0) ? (1) : (-1) )

// Type definitions

typedef ublas::vector<real> RealVectorL;
typedef ublas::matrix<real> RealMatrixL;
typedef ublas::vector<int> IntVectorL;
typedef ublas::matrix<int> IntMatrixL;


typedef ublas::bounded_vector<real, 3> RealVector;
typedef ublas::bounded_matrix<real, 3, 3> RealMatrix;
typedef ublas::bounded_vector<int, 3> IntVector;
typedef ublas::bounded_matrix<int, 3, 3> IntMatrix;

#endif /*SPH_GLOBAL_H_*/
