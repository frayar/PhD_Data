/**
 *
 * BIRCH++ 
 * Implements BIRCH Clustering
 *
 * class IRNG
 * Incremental construction of the Relative Neighbourhood Graph (RNG)
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-05-19
 * @version 1.0
 *
 */

#ifndef __RNG_H__
#define __RNG_H__

	#include "BIRCH.h"

	class RNG
	{
		// Methods
		public:
			static double Compute_RNG(double *pData, int nbData, int dimension, list<Edge> *edges, std::size_t *nb_edges);
			static double Compute_RNG(double *pData, int nbData, int dimension, std::map<__int64, double> *rng, std::size_t *nb_edges);
			static double IA_RNG_Insert(double *data, __int64 nbdata, __int64 dimension, std::map<__int64, double> *adjacencies, std::size_t *nb_edges);
			static double IA_RNG_Remove(double *data, __int64 nbdata, __int64 dimension, double *data_idx, std::map<__int64, double> *adjacencies, std::size_t *nb_edges);		//DTODO
			static double IA_RNG_Update(double *data, __int64 nbdata, __int64 dimension, double *data_idx, std::map<__int64, double> *adjacencies, std::size_t *nb_edges);		//DTODO
			static __int64 GetFarthest(std::map<__int64, double> neighbors);
			static __int64 GetNearest(std::map<__int64, double> neighbors);
			static double EuclidianDist( double *data, std::size_t dimension, std::size_t id1, std::size_t id2);

		// Attributes
		private:
			std::size_t				size;			/*! Number of objects in this cluster */
			std::vector<Object *>	elements;		/*! List of objects within this cluster */
	};

#endif