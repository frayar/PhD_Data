/**
 *
 * BIRCH++ 
 * Implements BIRCH Clustering
 *
 * main header
 * include, namespace, define, enum, struct, etc.
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-04-13
 * @version 1.0
 *
 */



// ----------------------------------------------------------------------------
//
// Memory leak detection
// https://msdn.microsoft.com/fr-fr/library/x98tx3cf%28v=vs.100%29.aspx
//
// ----------------------------------------------------------------------------
#ifdef _DEBUG
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
	#include <crtdbg.h>

	//#ifndef DBG_NEW
	//	#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
	//	#define new DBG_NEW
	//#endif
#endif  // _DEBUG
// ----------------------------------------------------------------------------

 #ifdef __linux__
 	#define __int64 int64_t
 #endif


#ifndef __BIRCH_H__
#define	__BIRCH_H__

	#pragma once

	#include <iostream>
	#include <fstream>
	#include <exception>
	#include <assert.h>
	#include <sys/time.h>
	//#include <ctime>
	#include <memory>

	#include <vector>
	#include <list>
	#include <map>
	#include<set>

	#include <math.h> 
	#include <limits>			// std::numeric_limits<float>::infinity()
	#include <algorithm>		// std::sort

	#include <omp.h>			// Use OMP

	using namespace std; 

	#define PAGE_SIZE					(4*1024) /* assuming 4K page */
	#define MAX_MEMORY					0
	#define THRESHOLD					20.0	/* 1.0 for NGA */ /* 25.0 for Wang */ /* 0.7 for MNIST*/ /* 20.0 for Castle */ /* 40.0 for NASA-APOD */
	#define INTERNAL_NODE_MAX_ENTRIES	50		/* 10 by default */
	#define LEAF_NODE_MAX_ENTRIES		50		/* 10 by default */	
	#define	NB_REPRESENTATIVES			7	
	#define IRNG_EPSILON				1		/* IRNG hypersphere magnifying value, between 0 and 1 */ 
	#define CURE						0		/* For farthest representative computation */
	#define iRNG_MINIMUM_OBJECTS		10000	/* If there is less that iRNG_MINIMUM_OBJECTS, it is faster to use parallelised RNG computation instead of incremental one */

	/*!
	 * \enum Label
	 * \brief Represents the position of the label if present.
	 */
	enum Label {
		START	= -1,	/**< Label at the first position */
		NONE	= 0,	/**< No label */
		END	    = 1		/**< Label at the last position */
	};


	/*!
	 * \struct Object [BIRCH.h]
	 * \brief Represent an object of the dataset
	 */
	struct Object {
		std::size_t id;				/**< ID of this data-point */
		std::string tree_id;		/**< ID of this data-point in the tree stucture */
		double* ptr;				/**< Pointer to object numeric signature */
		std::size_t dimension;		/**< Size of the signature */
		std::string label;			/**< Label of the object */
		std::string imagepath;		/**< Path of the corresponding image - Not assigned in the current version */ 
		std::size_t cluster_id;		/**< Id of the cluster (it may be changed as a std::string to reflect the hierachy, e.g. xx.xx.xx.xx */

		// Default constructor
		Object* clone()
		{
			Object *new_obj = new Object();
			new_obj->id = this->id;
			new_obj->tree_id = this->tree_id;
			new_obj->ptr = new double[dimension];
			std::copy(this->ptr, this->ptr + dimension, new_obj->ptr);
			new_obj->dimension = this->dimension;
			new_obj->label = this->label;
			new_obj->imagepath = this->imagepath;
			new_obj->cluster_id = this->cluster_id;

			return new_obj;
		}
	};

	
	/*!
	 * \struct Edge [BIRCH.h]
	 * \brief Represent an graph edge as a triplet
	 */
	struct Edge {
		std::size_t src;		/**< Source endpoint*/
		std::size_t tgt;		/**< Target endpoint */
		double weight;			/**< Weight */

		// Default constructor
		Edge() : 
			src(-1), tgt(-1), weight(0.0)  {}

		// Constructor
		Edge(std::size_t src_, std::size_t tgt_, double weight_) :
            src(src_), tgt(tgt_),weight(weight_) {}

		// Overload comparison method
		bool operator< (const Edge  &e) const 
		{
			return ( ( this->src < e.src) || (this->src == e.src && this->tgt < e.tgt) );
		}
	};

#endif
