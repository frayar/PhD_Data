/**
 * 
 * LayoutHRNG
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2016-08-17
 * @version 1.0
 *
 *
 * http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/graph_concepts.html
 * http://cpp.developpez.com/redaction/data/pages/users/gbdivers/boost-graph/
 * http://matthieu-brucher.developpez.com/tutoriels/cpp/boost/graph/implementation/
 *
 */


 #ifdef __linux__
 	#define __int64 int64_t
 #endif


#ifndef __LAYOUT_HRNG_H__
#define	__LAYOUT_HRNG_H__

	#pragma once

	#include <iostream>
	#include <fstream>
	#include <sstream>
	#include <exception>
	#include <assert.h>
	//#include <sys/time.h> For Fedora
	#include <ctime>
	#include <memory>
	#include <cstdio>
	#include <vector>
	#include <list>
	#include <map>
	#include<set>
	#include <math.h> 
	#include <limits>			// std::numeric_limits<float>::infinity()
	#include <algorithm>		// std::sort
	#include <omp.h>			// Use OMP

	// BOOST - JSON PARSING
	#include <boost/property_tree/ptree.hpp>
	#include <boost/property_tree/json_parser.hpp>

	// DOUBLE TO STRING CONVERSION
	#include <boost/lexical_cast.hpp>

	// BOOST - GRAPH STRUCTURE
	#include <boost/graph/graph_traits.hpp>
	#include <boost/graph/adjacency_list.hpp>
	#include <boost/graph/dijkstra_shortest_paths.hpp>
	#include <boost/property_map/property_map.hpp>
	#include <boost/property_map/transform_value_property_map.hpp>

	// NAMESPACE
	using namespace std; 
	using namespace boost::property_tree;

	struct VertexProperties 
	{ 
		std::string id;
	    double x;
	    double y;
		double sm_x;
		double sm_y;
	    VertexProperties() : id(""), x(0.0), y(0.0), sm_x(0.0), sm_y(0.0) {}
	    VertexProperties(std::string const& id_, double x_, double y_, double sm_x_, double sm_y_) : id(id_), x(x_), y(y_), sm_x(sm_x_), sm_y(sm_y_) {}
	};

	struct EdgeProperties 
	{ 
		std::string id;
		double weight;
	    EdgeProperties() : id(""), weight(0.0) {}
	    EdgeProperties(std::string const& id_, double weight_) : id(id_), weight(weight_) {}
	};

    typedef boost::adjacency_list<boost::vecS /* edges container */, boost::vecS /* vertices container */, boost::undirectedS , VertexProperties, EdgeProperties> Graph;
	typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef boost::graph_traits<Graph>::edge_descriptor Edge;

#endif
