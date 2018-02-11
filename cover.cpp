#include "cover.hpp"

CoverMax::CoverMax(vector<Curve>& curves, Kernel::FT& squared_radius)
{
	m_max_intersection			= 1;
	m_sqrd_radius			= squared_radius;
	m_distance_sqrd_radius	= Kernel::FT(DISTANCE_CALC * squared_radius);
	m_found					= false;

	insert(m_arr, curves.begin(), curves.end());
}

inline CGAL::Bounded_side bounded_side(const Kernel::Circle_2& circle, Kernel::FT sq_rad, const Point& p)
{
	Kernel::Point_2 center = circle.center();

	auto xd = Kernel::FT(center.x()) - p.x();
	auto yd = Kernel::FT(center.y()) - p.y();

	auto diff = xd*xd + yd*yd;

	switch (CGAL::compare(diff, sq_rad))
	{
	case CGAL::LARGER:
		return CGAL::ON_UNBOUNDED_SIDE;
	case CGAL::SMALLER:
		return CGAL::ON_BOUNDED_SIDE;
	//case CGAL::EQUAL:
	//	break;
	}
	return CGAL::ON_BOUNDARY;
}

void CoverMax::MapEdgesToCircles()
{
	unsigned int i = 1;

	for (Arr_with_hist_2::Curve_iterator cit = m_arr.curves_begin(); cit != m_arr.curves_end(); ++cit, i++)
	{
		DEBUG_PRINT("Circle number: " << i);

		//
		// going over all edges of a circle and creating maps
		//
		for (Arr_with_hist_2::Induced_edge_iterator ieit = m_arr.induced_edges_begin(cit); ieit != m_arr.induced_edges_end(cit); ++ieit)
		{
			if ((*ieit)->is_fictitious() == true)
			{
				DEBUG_PRINT("Fictitious edge: " << (*ieit)->source()->point() << "-->" << (*ieit)->target()->point());
			}
			else
			{
				map_edge_to_circle[*ieit] = *cit;
				map_edge_to_circle[(*ieit)->twin()] = *cit;
			}
		}
	}

	DEBUG_PRINT("");
}

void CoverMax::CreateCircleMap(Arr_with_hist_2::Curve_iterator& cit, circle_direction& cd)
{
	for (Arr_with_hist_2::Induced_edge_iterator ieit = m_arr.induced_edges_begin(cit); ieit != m_arr.induced_edges_end(cit); ++ieit)
	{
		if ((*ieit)->is_fictitious() == true)
		{
			DEBUG_PRINT("Fictitious edge: " << (*ieit)->source()->point() << "-->" << (*ieit)->target()->point());
		}

		else
		{
			Arr_with_hist_2::Vertex_const_handle curr_s, curr_t;
			Arr_with_hist_2::Halfedge_handle edge_in_wanted_direction;
			Arr_with_hist_2::X_monotone_curve_2 curve = (*ieit)->curve();

			if (((*ieit)->source()->point()) == curve.source())
			{
				curr_s = (*ieit)->source();
				curr_t = (*ieit)->target(); 
				edge_in_wanted_direction = *ieit;
			}
			else if (((*ieit)->target()->point()) == curve.source())
			{
				curr_s = (*ieit)->target();
				curr_t = (*ieit)->source();
				edge_in_wanted_direction = (*ieit)->twin();
			}

			cd[curr_s] = pair<Arr_with_hist_2::Vertex_const_handle, Arr_with_hist_2::Halfedge_handle>(curr_t, edge_in_wanted_direction);
		}
	}
}

void CoverMax::SetMaxVertex(Point& pnt)
{
	m_max_vertex = pnt;
	m_found = true;
}

void CoverMax::RunAlgorithm()
{
	unsigned int i = 1;

	//
	// creating a map between edges and circles curves
	//
	MapEdgesToCircles();

	//
	// iterating over the arrangement circles
	//
	for (Arr_with_hist_2::Curve_iterator cit = m_arr.curves_begin(); cit != m_arr.curves_end(); ++cit, i++)
	{
		DEBUG_PRINT("Circle number: " << i << endl << "----------------------");
		
		circle_direction map_vertex_linked_list;

		//
		// iterating over all edges of a circle, and
		// creating a map that produces a direction
		//
		CreateCircleMap(cit, map_vertex_linked_list);

		//
		// going over cit's vertices 2 times clockwise (by using the sorted vertices)
		//
		auto								first_key			= map_vertex_linked_list.begin()->first;
		Arr_with_hist_2::Halfedge_handle	current_edge;
		unsigned int						tangent_flag		= 0,
											max_intersection	= 0;
		auto								max_vertex			= first_key, 
											current_key			= first_key;
		map<ptrdiff_t, bool>				circles_appearance;
		
		for (int j = 0; j < CIRCLE_ITERATIONS; j++)
		{
			do
			{
				current_edge = map_vertex_linked_list[current_key].second;
				DEBUG_PRINT("ordered edge: " << current_edge->source()->point() << " --> " << current_edge->target()->point());

				if (current_key->degree() != POINT_ISNT_INTERSECTION)
				{
					DEBUG_PRINT("[+] intersection vertex detected");

					auto	first_incident_halfedge		= current_key->incident_halfedges(),
							current_incident_halfedge	= current_key->incident_halfedges();

					map<ptrdiff_t, bool> is_circle_already_seen;
					tangent_flag	= 0;
					vector<ptrdiff_t> circle_earse;

					//
					// going over all incident halfedges clockwise
					//
					do
					{
						// Note that the current halfedge is directed into the current_key
						Arr_with_hist_2::Halfedge_handle		current_incident_halfedge_handle	= static_cast<Arr_with_hist_2::Halfedge>(*current_incident_halfedge).twin()->twin();
						Arr_with_hist_2::Vertex_const_handle	incident_source_vertex				= current_incident_halfedge->source();

						Curve incident_cit = map_edge_to_circle[current_incident_halfedge_handle];
						DEBUG_PRINT("\n\t (" << incident_source_vertex->point() << ") on circle: " << incident_cit);

						//
						// check if we're entering, exiting or "touching" the incident circle,
						// according to the location of the next point on cit
						//
						Circle	cit_circ			= cit->supporting_circle(),
								incident_cit_circ	= incident_cit.supporting_circle();
						
						if (cit_circ == incident_cit_circ) 
						{
							DEBUG_PRINT("\t [-] incident edge is on current circle");
							continue;
						}
						if (is_circle_already_seen.find(incident_cit_circ.id()) != is_circle_already_seen.end())
						{
							DEBUG_PRINT("\t [-] circle was already spotted for incident edge");
							continue;
						}

						is_circle_already_seen[incident_cit_circ.id()] = true;
						
						DEBUG_PRINT("\t [+] incident edge is on different circle");
						DEBUG_PRINT("\t next point on cit: " << (map_vertex_linked_list[current_key].first)->point());
						CGAL::Bounded_side location = bounded_side(incident_cit_circ, m_sqrd_radius, (map_vertex_linked_list[current_key].first)->point());
						DEBUG_PRINT("\t event: " << location);

						// entering the incident circle
						if ((location == CGAL::ON_BOUNDED_SIDE) || (location == CGAL::ON_BOUNDARY)) 
						{
							circles_appearance[incident_cit_circ.id()] = true;
						}

						//
						// the circles are tangent OR exiting the incident circle
						// location == CGAL::ON_UNBOUNDED_SIDE
						//
						// check if the two circles are tangents by measuring the distance between their circle centers
						// and comparing its length to twice the size of the radius
						//
						else 
						{
							if (CGAL::squared_distance(incident_cit_circ.center(), cit_circ.center()) == m_distance_sqrd_radius)
							{
								// the circles are tangent
								tangent_flag = 1;
							}

							// exiting the incident circle
							else 
							{
								circle_earse.push_back(incident_cit_circ.id());
							}
						}
					} 
					while (++current_incident_halfedge != first_incident_halfedge);

					//
					// check whether we've found a new maximum, while
					// making sure to pick the lexicographic smallest vertex within a specific circle
					//
					
					unsigned int ac_size = static_cast<unsigned int>(circles_appearance.size() + tangent_flag);

					if (ac_size > max_intersection)
					{
						max_intersection = ac_size;
						max_vertex = current_key;
					}

					else if (ac_size == max_intersection)
					{
						Arr_with_hist_2::Traits_2::Compare_xy_2 comp;
						if (comp.operator()(current_key->point(), max_vertex->point()) == CGAL::SMALLER)
						{
							max_vertex = current_key;
						}
					}

					for (auto it = circle_earse.begin(); it != circle_earse.end(); ++it)
					{
						if (circles_appearance.find(*it) != circles_appearance.end())
						{
							circles_appearance.erase(*it);
						}
					}
				}
				else
				{
					DEBUG_PRINT("[-] none intersection vertex detected");

					if (!m_found)
					{
						auto pnt = current_key->point();
						Arr_with_hist_2::Traits_2::Compare_xy_2 comp;
						if (comp.operator()(pnt, m_max_vertex) == CGAL::SMALLER)
						{
							m_max_vertex = pnt;
						}
					}
				}

				current_key = map_vertex_linked_list[current_key].first;
			} 
			while (first_key != current_key);

			DEBUG_PRINT("");
		}

		//
		// update global max intersection and vertex across all circles
		//
		max_intersection++;
		auto max_point = max_vertex->point();
		
		if (max_intersection > m_max_intersection)
		{
			m_max_intersection = max_intersection;
			SetMaxVertex(max_point);
		}

		else if (max_intersection == m_max_intersection)
		{
			Arr_with_hist_2::Traits_2::Compare_xy_2 comp;
			if (comp.operator()(max_point, m_max_vertex) == CGAL::SMALLER)
			{
				SetMaxVertex(max_point);
			}
		}
	}
}

