/**
 *
 * BIRCH++ 
 * Implements BIRCH Clustering
 *
 * class CFNode
 * Node structure of an CFTree
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-04-13
 * @version 1.0
 *
 */

#ifndef __CF_NODE_H__
#define	__CF_NODE_H__

	#include "BIRCH.h"

	// Forward declaration
	class CFEntry;

	class CFNode 
	{
		// Methods
		public:
			// Constructor
			CFNode();

			// Methods
			void Add( CFEntry& new_entry );
			void Insert( CFEntry& new_entry, std::size_t position );
			std::size_t Replace( CFEntry& old_entry, CFEntry& new_entry );
			std::size_t GetHeight();
			std::size_t GetObjectsNumber();
			bool IsFull();
			bool IsEmpty();
			void AssignInitialID();
			void ReassignID();
			void CreateHierarchicalRNG();
			void UpdateRNG();

			// Getter/Setter
			std::size_t GetSize();
			CFEntry* GetEntries();
			std::string GetID();
			void SetID(std::string id_);
			bool GetLeafAttribute();
			void SetLeafAttribute(bool is_leaf_);
			CFNode* GetPrevious();
			void SetPrevious(CFNode* prev_);
			CFNode* GetNext();
			void SetNext(CFNode* next_);
			CFNode* GetFirstLeaf();
			void SetFirstLeaf(CFNode* first);
			CFNode* GetLastLeaf();
			void SetLastLeaf(CFNode* last);
			std::list<Edge>* GetRNG();
			std::size_t GetNbEdges();

			
			// Notification
			std::size_t Print(std::size_t nb_tab);


		// Attributes
		private:
			std::string					id;							/*! ID of the CFNode */
			std::size_t					size;						/*! Number of CFEntries contained in this CFNode */
			CFEntry*					entries;					/*! Pointer to an array of this CFNode entries */
			bool						is_leaf;					/*! Boolean that indicates if the node is a leaf */
			CFNode*						prev;						/*! Previous CFNode */
			CFNode*						next;						/*! Next CFNode */
			CFNode*						first_leaf;					/*! First leaf in the subtree of this node */
			CFNode*						last_leaf;					/*! Last leaf in the subtree of this node */
			std::list<Edge>*			rng;						/*! RNG of the entries of this node */
			std::size_t					nb_edges;					/*! Number of edge in the rng */
	};

#endif
