/**
 *
 * BIRCH++ 
 * Implements BIRCH Clustering
 *
 * class CFEntry
 * Summarize a cluster of data-points
 * According to BIRCH paper, an CFEntry is given by a triple (N, LS, SS)
 * - N is the number of data-points in the cluster
 * - LS is the linear sum of the N data-points (it is a vector)
 * - SS is the square sum of the N data-points (it is a value)
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-04-13
 * @version 1.0
 *
 */

#ifndef __CF_ENTRY_H__
#define	__CF_ENTRY_H__

	#include "BIRCH.h"
	#include "Utils.h"
	#include "CFNode.h"
	#include "Cluster.h"

	class CFEntry
	{
		// Methods
		public:
			// Constructors
			CFEntry();
			CFEntry(Object* dataPoint);
			CFEntry(CFNode* node);

			// Operators overload
			void operator+=( const CFEntry& right_side_entry );
			CFEntry operator+( const CFEntry& right_side_entry );
			void operator-=( const CFEntry& right_side_entry );

			// Getter/Setter
			std::string GetID();
			void SetID(std::string id_);
			CFNode* GetChild();
			Object* GetObject();
			std::size_t GetObjectsNumber();
			std::size_t GetDimension();
			static void SetDimension(std::size_t dimension_);
			Cluster* GetClusters();
			std::vector<Object *> GetFarRepresentatives();
			std::vector<Object *> GetNearRepresentatives();

			// Others methods
			void GetCentroid(double *centroid);
			bool HasChild();
			bool IsInLeaf();
			void ComputeRepresentatives();
			void UpdateRepresentatives();

			// Notification
			std::size_t Print(std::size_t nb_tab);

			// Cluster attributes
			static Object Centroid( const CFEntry& entry );
			static double Radius( const CFEntry& entry );
			static double Diameter( const CFEntry& entry );

			// Distances
			static double DistD0( const CFEntry& left_side_entry, const CFEntry& right_side_entry );
			static double DistD1( const CFEntry& left_side_entry, const CFEntry& right_side_entry );
			static double DistD2( const CFEntry& left_side_entry, const CFEntry& right_side_entry );
			static double DistD3( const CFEntry& left_side_entry, const CFEntry& right_side_entry );
			static double DistD4( const CFEntry& left_side_entry, const CFEntry& right_side_entry );
			static double EuclidianDist( const std::size_t dimension, const double* left_side_data_point, const double* right_side_data_point );


		// Attributes
		private:
			std::string					id;						/*! ID of the CFEntry */
			Object*						object;					/*! Pointer to the related Object */
			static std::size_t			dim;					/*! Dimension of a data-point */
			std::size_t					n;						/*! Number of data-points */
			double*						lin_sum;				/*! Linear sum of each dimension of n data-points */
			double						sq_sum;					/*! Square sum of n data-points */
			CFNode*						child;					/*! Pointer to a child node, null if the CFEntry belongs to a leaf */ 
			//shared_ptr<CFNode>		parent;					/*! Pointer to the node containing this CFEntry */
			Cluster*					clusters;				/*! Pointer to list of clusters within this entry (one by default if is leaf) */
			std::vector<Object *>		near_representatives;	/*! Representatives of the current node subtree : nearest neighbors of the medoid */
			std::vector<Object *>		far_representatives;	/*! Representatives of the current node subtree : farthest neighbors of the medoid (CURE-like) */

	};

#endif
