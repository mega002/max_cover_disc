#include <boost/timer.hpp>
#include "cover.hpp"

//
// globals
//
const char *usage = "USAGE: <squared radius> <input file>";

//
// heplers
//
void execute(curves& curves, Kernel::FT& radius, unsigned int& cover, Point& point);

std::ostream& operator<<(std::ostream& os, const CoordNT& c);

std::ostream& operator<<(std::ostream& os, const Point& p);

//
// main
//
int main(int argc, char* argv[])
{
	Kernel::FT squared_radius;
	Point vertex;
	unsigned int cover = 0;

	if (argc < 2) {
		std::cerr << "The radius is missing!" << std::endl;
		return 1;
	}
	std::istringstream r_square_stream(argv[1]);
	r_square_stream >> squared_radius;

	const char* filename = (argc > 2) ? argv[2] : "points.txt";
	std::ifstream is;
	is.open(filename);
	if (!is.is_open()) {
		std::cerr << "Failed to open " << filename << "!" << std::endl;
		return 1;
	}
	size_t n;
	is >> n;
	std::vector<Curve> curves(n);
	for (auto i = 0; i < n; ++i) {
		Kernel::FT x, y;
		is >> x >> y;
		auto p = Rational_point(x, y);
		curves[i] = Curve(Circle(p, squared_radius));
	}

	boost::timer timer;

	execute(curves, squared_radius, cover, vertex);

	DEBUG_PRINT("comparison: " << vertex.equals(Traits::Point_2(0, 0)));

	cout << "Execution time: " << timer.elapsed() << endl;
	cout << "The maximum covering disc covers " << cover << " points" << endl;
	cout << "The maximum covering disc is centered at: " << vertex << endl;
	
	return ERR_SUCCESS;
}

inline void execute(curves& curves, Kernel::FT& radius, unsigned int& cover, Point& point)
{
	// generate covermax instance for calculations
	CoverMax cm = CoverMax(curves, radius);

	// execute the algorithm
	cm.RunAlgorithm();

	// return the results
	cover = cm.GetMaxCounter();
	point = cm.GetMaxVertex();
}

std::ostream& operator<<(std::ostream& os, const CoordNT& c)
{
	CGAL::Rational_traits<Number_type> rt;

	std::cout << "("
		<< rt.numerator(c.a0()) << "/" << rt.denominator(c.a0()) << ","
		<< rt.numerator(c.a1()) << "/" << rt.denominator(c.a1()) << ","
		<< rt.numerator(c.root()) << "/" << rt.denominator(c.root()) << ")";
	return os;
}

std::ostream& operator<<(std::ostream& os, const Point& p)
{
	std::cout << "(" << p.x() << "," << p.y() << ")";
	return os;
}