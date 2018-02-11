#pragma once

#include "globals.hpp"

class CoverMax
{
public:
	CoverMax(vector<Curve>& curves, Kernel::FT& squared_radius);

	//
	// creates a map between vertices and its circles
	//
	void RunAlgorithm();

	//
	// returns the maximum intersection found by the algorithm
	//
	unsigned int GetMaxCounter() { return m_max_intersection; }

	//
	// returns a center point based on vertex with maximum intersections
	//
	Point& GetMaxVertex() { return m_max_vertex; }

	void SetMaxVertex(Point& pnt);

protected:
	void MapEdgesToCircles();
	void CreateCircleMap(Arr_with_hist_2::Curve_iterator& cit, circle_direction& cd);

	//
	// members
	//
	Kernel::FT m_sqrd_radius;
	Kernel::FT m_distance_sqrd_radius;
	Arr_with_hist_2 m_arr;
	Point m_max_vertex;
	unsigned int m_max_intersection;
	bool m_found;
	map<Arr_with_hist_2::Halfedge_handle, Curve> map_edge_to_circle;
};
