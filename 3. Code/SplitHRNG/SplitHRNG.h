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

	// NAMESPACE
	using namespace std; 
	using namespace boost::property_tree;


#endif
