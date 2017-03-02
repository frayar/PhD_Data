/**
 *
 * BIRCH++ 
 * Implements BIRCH Clustering
 *
 * class Cluster
 * Actual cluster, i.e. group of homegeneous data-points
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-04-13
 * @version 1.0
 *
 */

#ifndef __CLUSTER_H__
#define __CLUSTER_H__

	#include "BIRCH.h"

	class Cluster
	{
		// Methods
		public:
			Cluster();
			void operator+=(Cluster* right_side_cluster);
			void AddElement(Object* element);
			void AddElements(std::vector<Object *> elements_);
			std::size_t GetSize();
			std::vector<Object *> GetElements();
			void CreateRNG(std::string entry_id);
			void UpdateRNG(Cluster* new_cluster);
			void ClearRNG();
			std::string GetID();
			void SetID(std::string id_);
			std::map<__int64, double>* GetRNG();
			std::size_t GetNbEdges();

		// Attributes
		private:
			std::string					id;				/*! ID of the cluster in the Tree. Same as its parent CFEntry */
			std::size_t					size;			/*! Number of objects in this cluster */
			std::vector<Object *>		elements;		/*! List of objects within this cluster */
			std::map<__int64, double>*	rng;			/*! RNG of the element of this cluster */
			std::size_t					nb_edges;		/*! Number of edge in the rng */
	};

#endif
