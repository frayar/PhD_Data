/**
 *
 * BIRCH++ 
 * Implements BIRCH Clustering
 *
 * class CFTree ( clustering feature tree )
 * B-tree-like data structure consisting of summarized clusters
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-04-13
 * @version 1.0
 *
 */
#ifndef __CF_TREE_H__
#define	__CF_TREE_H__

	#include "BIRCH.h"
	#include "CFNode.h"
	#include "CFEntry.h"

	class CFTree
	{
		// Methods
		public:
			// Constructors
			CFTree(std::size_t dimension);
			CFTree(std::size_t dimension, double threshold_, std::size_t internal_node_max_entries_, std::size_t leaf_max_entries_);

			// Getter/Setter
			CFNode* GetRoot();
			static std::size_t GetInternalNodeMaxEntries();
			static std::size_t GetLeafNodeMaxEntries();
			std::size_t GetHeight();
			std::vector<std::vector<CFEntry*>>* GetRngsNodes();
			std::vector<std::vector<Edge>>* GetRngsEdges();
			CFNode* GetDummyFirstLeaf();

			// Tree manipulation methods
			void Insert(CFEntry& entry);
			void Insert(CFNode* node, CFEntry& entry, bool &bsplit);
			void SplitRoot( CFEntry& entry);
			void SplitNode( CFNode& node, CFEntry& closest_entry, CFEntry& entry, bool& bsplit );
			void RedistributeEntries( vector<CFEntry> *entries, std::pair<CFEntry*,CFEntry*> farthest_pair, CFEntry& new_entry_left, CFEntry& new_entry_right );
			void UpdatePointers( CFNode* new_node_left, CFNode* new_node_right );
			void Rebuild( bool extend = true );
			std::size_t ComputeHeight();
			std::vector<Cluster *> GetClusters();
			std::vector<Cluster *> GetPrunedClusters();
			void AssignInitialID();
			void ComputeRepresentatives();
			void CreateHierarchicalRNG();
			void CreateMultilevelRNG();
			void InsertAndUpdate(CFEntry& entry);
			void InsertAndUpdate(CFNode* node, CFEntry& entry, bool &bsplit);
			void SplitRootAndUpdate( CFEntry& entry);
			void SplitNodeAndUpdate( CFNode& node, CFEntry& closest_entry, CFEntry& entry, bool& bsplit );

			// Distance related methods
			double AverageDistanceOfClosestPairLeafEntries();
			CFEntry* FindClosestEntry( CFEntry* entry, CFNode* node );
			void FindFarthestPair( std::vector<CFEntry> *entries, std::pair<CFEntry*, CFEntry*>& farthest_pair);

			// Notification
			void PrintLeaves();
			std::size_t Print();


		// Attributes
		private:
			CFNode*									root;						/*! Pointer to rooot */
			CFNode*									dummy_first_leaf;			/*! Pointer to first leaf */

			static std::size_t						dim;						/*! Dimension of the data */
			static std::size_t						internal_node_max_entries;	/*! Internal node max entries B */
			static std::size_t						leaf_max_entries;			/*! Leaf max entries L */
			double									threshold;					/*! Distance threshold T */

			std::size_t								node_cnt;					/*! Memory usage (for monitoring) */
			std::size_t								mem_limit;					/*! Max memory allowed */
			std::size_t								height;						/*! Height of the tree */
			std::vector<std::vector<Edge>>			rngs_edges;					/*! Vector of multilevel rng */
			std::vector<std::vector<CFEntry*>>		rngs_nodes;					/*! Vector of multilevel rng nodes. An RNG node corresponds to the CFEntries */
			std::vector<std::size_t>				rngs_nb_nodes;				/*! Vector of multilevel rng nodes number*/
			std::vector<std::size_t>				rngs_nb_edges;				/*! Vector of multilevel rng edges number*/

	};

#endif
