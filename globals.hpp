#pragma once

#include <fstream>
#include <sstream>
#include <ostream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_with_history_2.h>

using namespace std;

//
// CGAL types
//

typedef CGAL::Exact_predicates_exact_constructions_kernel	Kernel;
typedef Kernel::FT											Number_type;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>			Traits;
typedef Traits::CoordNT										CoordNT;
typedef Traits::Point_2										Point;
typedef Traits::Curve_2										Curve;
typedef Traits::Rational_point_2							Rational_point;
typedef Traits::Rational_segment_2							Segment;
typedef Traits::Rational_circle_2							Circle;
typedef CGAL::Arrangement_with_history_2<Traits>			Arr_with_hist_2;

// 
// CoverMax types
//
typedef Point																					vertex;
typedef vector<Curve>																			curves;
typedef map<Arr_with_hist_2::Vertex_const_handle, pair<Arr_with_hist_2::Vertex_const_handle, 
			Arr_with_hist_2::Halfedge_handle>>													circle_direction;

#undef _INTERNAL_DEBUG
#ifdef _INTERNAL_DEBUG
#define DEBUG_PRINT(x)	{ cout << x << endl; }
#else
#define DEBUG_PRINT(x)
#endif // _DEBUG


const unsigned int	POINT_ISNT_INTERSECTION = 2;
const int			ERR_SUCCESS				= 0;
const int			ERR_FAILURE				= 1;

const int			ARGC_VALID_COUNT		= 3;
const int			ARGC_FILE_POS			= 2;
const int			ARGC_RADIUS_POS			= 1;

const double		DISTANCE_CALC			= 4;
const unsigned int	SINGLE_CIRCLE			= 1;
const unsigned int	CIRCLE_ITERATIONS		= 2;