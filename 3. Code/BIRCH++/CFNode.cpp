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

#include "Utils.h"
#include "RNG.h"
#include "CFNode.h"
#include "CFEntry.h"
#include "CFTree.h"
#include "BIRCH.h"

/*****************************************************************************/
/*  CONSTRUCTORS                                                             */
/*****************************************************************************/

/*!
 *	\brief Constructor (default)
 */
CFNode::CFNode()
{
	// Initialize attributes
	size = 0;
	is_leaf = false;
	prev = NULL;
	next = NULL;
	first_leaf = this;
	last_leaf = this;
	nb_edges = 0;

	// Create an array of CFEntry - N.B. A CfNode should fit in a memory page 
	//entries = new CFEntry[(PAGE_SIZE - ( sizeof(CFNode*)*2 /* 2 leaf node pointers */ + sizeof(std::size_t) /* size */ )) / sizeof(CFEntry)/*max_entries*/];	
	entries = new CFEntry[INTERNAL_NODE_MAX_ENTRIES];
}


/*****************************************************************************/
/*  GETTER/SETTER                                                            */
/*****************************************************************************/

/*!
 *	\brief Get node size
 *
 *	\return Node size
 */
std::size_t CFNode::GetSize()
{
	return size;
}


/*!
 *	\brief Get pointer to entries array
 *
 *	\return Pointer to entries array
 */
CFEntry* CFNode::GetEntries()
{
	return entries;
}


/*!
 *	\brief Indicate if the node is a leaf 
 *
 *	\return True id the node is a leaf, else False
 */
std::string CFNode::GetID()
{
	return id;
}


/*!
 *	\brief Define if the node is a leaf 
 *
 *	\param is_leaf_ : leaf boolean status
 *
 */
void CFNode::SetID(std::string id_)
{
	id =  id_;
}


/*!
 *	\brief Indicate if the node is a leaf 
 *
 *	\return True id the node is a leaf, else False
 */
bool CFNode::GetLeafAttribute()
{
	return is_leaf;
}


/*!
 *	\brief Define if the node is a leaf 
 *
 *	\param is_leaf_ : leaf boolean status
 *
 */
void CFNode::SetLeafAttribute(bool is_leaf_)
{
	is_leaf = is_leaf_;
}


/*!
 *	\brief Get a pointer to the previous node 
 *
 *	\return Pointer to the previous node
 */
CFNode* CFNode::GetPrevious()
{
	return prev;
}


/*!
 *	\brief Set the previous node pointer
 *
 *	\param	prev_ : pointer to the previous node
 */
void CFNode::SetPrevious(CFNode* prev_)
{
	prev = prev_;
}


/*!
 *	\brief Get a pointer to the next node 
 *
 *	\return Pointer to the next node
 */
CFNode* CFNode::GetNext()
{
	return next;
}


/*!
 *	\brief Set the next node pointer
 *
 *	\param	prev_ : pointer to the next node
 */
void CFNode::SetNext(CFNode* next_)
{
	next = next_;
}


/*!
 *	\brief Get a pointer to the first leaf of the node subtree
 *
 *	\return Pointer to the first leaf of the node subtree
 */
CFNode* CFNode::GetFirstLeaf()
{
	return first_leaf;
}


/*!
 *	\brief Set the pointer of the first leaf of the node subtree
 *
 *	\param	prev_ : pointer to thefirst leaf of the node subtree
 */
void CFNode::SetFirstLeaf(CFNode* first)
{
	first_leaf = first;
}


/*!
 *	\brief Get a pointer to the first leaf of the node subtree
 *
 *	\return Pointer to the first leaf of the node subtree
 */
CFNode* CFNode::GetLastLeaf()
{
	return last_leaf;
}


/*!
 *	\brief Set the pointer of the first leaf of the node subtree
 *
 *	\param	prev_ : pointer to thefirst leaf of the node subtree
 */
void CFNode::SetLastLeaf(CFNode* last)
{
	last_leaf = last;
}

/*!
 *	\brief Get a pointer to the first leaf of the node subtree
 *
 *	\return Pointer to the first leaf of the node subtree
 */
std::list<Edge>* CFNode::GetRNG()
{
	return rng;
}

/*!
 *	\brief Get a pointer to the first leaf of the node subtree
 *
 *	\return Pointer to the first leaf of the node subtree
 */
std::size_t CFNode::GetNbEdges()
{
	return nb_edges;
}

			



/*****************************************************************************/
/*  NODE MANIPULATION METHODS                                                */
/*****************************************************************************/

/*!
 *	\brief Add new CFEntry to this CFNode
 *
 *	\param newEntry	: new entry
 *
 */
void CFNode::Add( CFEntry& new_entry )
{
	assert( size < (GetLeafAttribute() ? CFTree::GetLeafNodeMaxEntries() : CFTree::GetInternalNodeMaxEntries()) );
	entries[size++] = new_entry;
}


/*!
 *	\brief Insert new CFEntry to this CFNode
 *
 *	\param newEntry	: new entry
 *	\param pos : position of the new entry
 *
 */
void CFNode::Insert( CFEntry& new_entry, std::size_t position )
{
	assert( size < (GetLeafAttribute() ? CFTree::GetLeafNodeMaxEntries() : CFTree::GetInternalNodeMaxEntries()) );

	// Insert the new entry at the correct position
	for (std::size_t i = size; i > position; i--)
		entries[i] = entries[i-1];

	// Insert the entry at the right position
	entries[position] = new_entry;

	// Increment number of entries
	size++;
}

/*!
 *	\brief Replace old CFEntry as new CFEntry
 *
 *	\param oldEntry : old entry to replace
 *	\param newEntry	: new entry
 *
 */
std::size_t CFNode::Replace( CFEntry& old_entry, CFEntry& new_entry )
{
	for( std::size_t i = 0 ; i < size ; i++ )
	{
		if( &entries[i] == &old_entry )
		{
			entries[i] = new_entry;
			return i;
		}
	}
	// should never come here
	assert(false);
}


/*!
 *	\brief Get height of the CFNode in the CTFree 
 *
	\return Node height
 */
std::size_t CFNode::GetHeight()
{
	std::size_t max_height = 0;
	// If leaf
	if (GetLeafAttribute())
		return 0;
	// If not leaf, go deeper in the subtree
	else
		for (std::size_t i=0; i<size; i++)
			max_height = std::max( max_height, (entries[i].GetChild())->GetHeight() );

	return ( 1 + max_height);
}

/*!
 *	\brief Get number of items within the leaf node entries
 *
 *	\return Number of objects within this leaf entries
 */
std::size_t CFNode::GetObjectsNumber()
{
	assert(is_leaf);

	std::size_t nb = 0;
	for (std::size_t i = 0; i<size; i++)
		nb += entries[i].GetObjectsNumber();

	return nb;
}

/*!
 *	\brief CFNode is full, no more CFEntries can be in 
 */
bool CFNode::IsFull()
{
	return size == (GetLeafAttribute() ? CFTree::GetLeafNodeMaxEntries() : CFTree::GetInternalNodeMaxEntries());
}

/*!
 *	\brief CFNode has nothing
 */
bool CFNode::IsEmpty()
{
	return size == 0;
}


/*****************************************************************************/
/*  RNG                                                                      */
/*****************************************************************************/

/*!
 *	\brief Batch mode: travel along the tree, assign initial IDs, 
 *
 */
void CFNode::AssignInitialID()
{
	// Assign IDs of node entries
	for (std::size_t i=0; i<size; i++)
	{
		std::string tmpID = id; tmpID += "."; tmpID += utils_to_string(i);
		entries[i].SetID(tmpID);
	}

	// Go to next level
	for (std::size_t i=0; i<size; i++)
		if (!is_leaf) 
			entries[i].GetChild()->AssignInitialID();
}

/*!
 *	\brief Dynamic mode: travel along the tree, assign initial IDs, 
 *
 */
void CFNode::ReassignID()
{
	// Assign IDs of node entries
	for (std::size_t i=0; i<size; i++)
		entries[i].SetID(id + "." + utils_to_string(i));

	// Go to next level
	for (std::size_t i=0; i<size; i++)
		if (!is_leaf) 
			entries[i].GetChild()->ReassignID();
}


/*!
 *	\brief Batch mode :
 *		   (i) create RNG of its entries
 *         (ii) process next level
 *
 */
void CFNode::CreateHierarchicalRNG()
{
	// Create RNG of node entries
	UpdateRNG();

	// Go to next level
	for (std::size_t i=0; i<size; i++)
		if (!is_leaf) 
			entries[i].GetChild()->CreateHierarchicalRNG();
		else // We shall create RNG at cluster level
			entries[i].GetClusters()->CreateRNG(entries[i].GetID());
}


/*!
 *	\brief Create/Update RNG of its entries centroid
 *
 */
void CFNode::UpdateRNG()
{
	// Proceed only if at least two child entries
	if (size > 1)
	{
		// Reset list if it exists
		if (nb_edges > 0)
		{
			for (std::size_t i=0; i<size-1; i++)
				rng[i].clear();
			nb_edges = 0;
			delete[] rng;
		}
	
		//Initialize edge list structure
		rng = new list<Edge>[size];

		// Get dimension of the data
		std::size_t dimension = entries[0].GetDimension();

		// Create an array of pointers to centroids
		double *pData = new double[size*dimension];
		double *ptr = pData;

		for (std::size_t i=0; i<size; i++)
		{
			ptr = pData + i*dimension;
			entries[i].GetCentroid(ptr);
		}

		// Call RNG.cpp function
		double time = RNG::Compute_RNG(pData, size, dimension, rng, &nb_edges);

		// Debug - 04/11/2015
		//if (nb_edges == 0)
		//	bool stop = true;

		// Free memory
		delete[] pData;
	}
}




/*****************************************************************************/
/*  NOTIFICATION                                                             */
/*****************************************************************************/

/*!
 *	\brief Print a CFNode structure in the console
 */
std::size_t CFNode::Print(std::size_t nb_tab)
{
	std::size_t cnt = 0;
	std::string s_leaf = is_leaf? " leaf=true" : "";

	std::string tabulation = std::string(nb_tab, '\t');

	std::cout << tabulation.c_str() << "<node id=" << id.c_str() << s_leaf.c_str() << " >" << std::endl;

	for (std::size_t i=0; i<size; i++)
	{
		cnt += entries[i].Print(nb_tab+1);
	}

	std::cout << tabulation.c_str() << "</node>" << std::endl;

	return cnt;

}
