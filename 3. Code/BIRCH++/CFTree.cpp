/**
 *
 * BIRCH++ 
 * Implements BIRCH Clustering
 *
 * class CFTree ( clustering feature tree )
 * Actual cluster, i.e. group of homegeneous data-points
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-04-13
 * @version 1.0
 *
 */

#include "Utils.h"
#include "RNG.h"
#include "CFTree.h"

// Static attributes initialisation
std::size_t CFTree::dim = 0;
std::size_t CFTree::internal_node_max_entries = INTERNAL_NODE_MAX_ENTRIES;
std::size_t CFTree::leaf_max_entries = LEAF_NODE_MAX_ENTRIES;


/*****************************************************************************/
/*  CONSTRUCTORS                                                             */
/*****************************************************************************/

/*!
 *
 *	\brief Constructor (default)
 *		   Initially, root is considered as the first leaf
 *
 *	\param dimension : dimension of the data
 *
 */
CFTree::CFTree(std::size_t dimension)
{
	CFEntry::SetDimension(dimension);
	root = new CFNode();
	root->SetID("0");
	root->SetLeafAttribute(true);
	dummy_first_leaf = root;
	node_cnt = 1;
	dim = dimension;
	threshold = THRESHOLD;
	mem_limit = MAX_MEMORY;

	if ( CFTree::internal_node_max_entries == 0 || CFTree::leaf_max_entries == 0)
	{
		CFTree::internal_node_max_entries = ( PAGE_SIZE - ( sizeof(CFNode*)*2 + sizeof(std::size_t) ) ) / sizeof(CFEntry);
		CFTree::leaf_max_entries = ( PAGE_SIZE - ( sizeof(CFNode*)*2 + sizeof(std::size_t) ) ) / sizeof(CFEntry);
	}
}


/*!
 *	\brief Constructor
 *		   Initially, root is considered as the first leaf
 *
 *	\param dimension : dimension of the data
 *	\param threshold_					: distance threshold T
 *	\param internal_node_max_entries_	: internal node max entries B
 *	\param leaf_max_entries_			: leaf max entries L
 *
 */
CFTree::CFTree(std::size_t dimension, double threshold_, std::size_t internal_node_max_entries_, std::size_t leaf_max_entries_)
{

	CFEntry::SetDimension(dimension);
	root = new CFNode();
	root->SetID("0");
	root->SetLeafAttribute(true);
	dummy_first_leaf = root;
	node_cnt = 1;
	dim = dimension;
	threshold = threshold_;
	mem_limit = MAX_MEMORY;
	CFTree::internal_node_max_entries = internal_node_max_entries_;
	CFTree::leaf_max_entries = leaf_max_entries_;

}

/*****************************************************************************/
/*  GETTER/SETTER                                                            */
/*****************************************************************************/

/*!
 *	\brief Get the maximum allowed entries for an leaf node
 *
 *	\return Maximum allowed entries for an leaf node
 */
CFNode* CFTree::GetRoot()
{
	return root;
}

/*!
 *	\brief Get the maximum allowed entries for an internal node
 *
 *	\return Maximum allowed entries for an internal node
 */
std::size_t CFTree::GetInternalNodeMaxEntries()
{
	return internal_node_max_entries;
}

/*!
 *	\brief Get the maximum allowed entries for an leaf node
 *
 *	\return Maximum allowed entries for an leaf node
 */
std::size_t CFTree::GetLeafNodeMaxEntries()
{
	return leaf_max_entries;
}

/*!
 *	\brief Get the height of the tree
 *
 *	\return Height of the tree
 */
std::size_t CFTree::GetHeight()
{
	return height;
}

/*!
 *	\brief Get a pointer to the list of RNGs nodes
 *
 *	\return Pointer to the list of RNGs nodes
 */
std::vector<std::vector<CFEntry*>>* CFTree::GetRngsNodes()
{
	return &rngs_nodes;
}

/*!
 *	\brief Get a pointer to the list of RNGs edges
 *
 *	\return Pointer to the list of RNGs edges
 */
std::vector<std::vector<Edge>>* CFTree::GetRngsEdges()
{
	return &rngs_edges;
}

/*!
 *	\brief Get a pointer to the dummy first leaf
 *
 *	\return Pointer to the dummu first leaf
 */
CFNode* CFTree::GetDummyFirstLeaf()
{
	return dummy_first_leaf;
}

/*****************************************************************************/
/*  TREE MANIPULATION METHODS                                                */
/*****************************************************************************/

/*!
 *	\brief Insert a CFEntry in the CFTree
 *
 *	\param entry : entry to insert in the CFTree
 *
 */
void CFTree::Insert(CFEntry& entry)
{
	bool bsplit;

	// Insert new entry in the root
	Insert(root, entry, bsplit);

	// If a split has been backpropagated up to the root, and this last one need to be split
	if( bsplit )
		SplitRoot( entry );

	// OutOfMemory management
	std::size_t curr_mem = node_cnt * sizeof(CFNode);
	if( mem_limit > 0 && node_cnt * sizeof(CFNode) > mem_limit )
	{
		Rebuild();
	}

	return;
}


/*!
 *	\brief Insert a CFEntry in the given CFNode
 *
 *	\param node  : node to consider
 *	\param entry : entry to insert in the CFTree
 *
 */
void CFTree::Insert(CFNode* node, CFEntry& entry, bool &bsplit)
{
	// Empty node (for first insertion in the tree)
	if (node->IsEmpty())
	{
		// Add the entry to the node
		node->Add(entry);
		// No need to split anything
		bsplit = false;
		// Return with success
		return;
	}

	// Find the closest entry among node children
	CFEntry *closest_entry = FindClosestEntry( &entry, node );

	// TODO: Here, we could add a contribution
	// If distance between the new data-point and it closest entry is too large w.r.t. a given threslhold
	// We'll rather create a new entry at this node level, instead of going deeper in the tree
	// This may avoid split operations

	// If internal node
	if (closest_entry->HasChild())
	{
		// Recursively insert the entry in the child node
		Insert( closest_entry->GetChild(), entry, bsplit );

		// Update CFEntry if no split
		if( !bsplit )
			*closest_entry += entry;
		// else split here
		else
			SplitNode( *node, *closest_entry, entry, bsplit );
	}
	// Else leaf
	else {
		// if the leaf node can absorb the new entry
		if ( CFEntry::DistD0(*closest_entry, entry) < threshold)
		{
			// Insert entry and update the leaf in the same time
			*closest_entry += entry;
			// No need to split, backpropagate this information to update at parent-level
			bsplit = false;
		}
		// else add new_entry if possible
		else if( node->GetSize() < CFTree::leaf_max_entries )
		{
			node->Add(entry);
			bsplit = false;
		}
		// if not possible, notify the split situation to the parent who will handle the actual split operation
		else	
		{
			bsplit = true;
		}
	}

	return;
}


/*!
 *	\brief Split the root
 *
 *	\param entry : new entry to insert while splitting
 *
 */
void CFTree::SplitRoot( CFEntry& entry )
{
	// Get the list of entries to consider (i.e. root entries plus the new entry)
	std::vector<CFEntry> entries; 
	entries.reserve( CFTree::internal_node_max_entries + 1 );
	for( std::size_t i = 0 ; i < CFTree::internal_node_max_entries ; i++ )
		entries.push_back((root->GetEntries())[i]);
	entries.push_back(entry);

	// Find the farthest entry pair
	std::pair<CFEntry*, CFEntry*> farthest_pair;
	FindFarthestPair( &entries, farthest_pair );

	// Create two new nodes
	CFNode* new_node_left = new CFNode();
	CFNode* new_node_right = new CFNode();

	// Create two entries from the two new nodes
	CFEntry new_entry_left( new_node_left );
	CFEntry new_entry_right( new_node_right );

	// Redistribute old entries to new entries
	RedistributeEntries(&entries, farthest_pair, new_entry_left, new_entry_right);

	// Get leaf status of the root
	bool root_is_leaf = root->GetLeafAttribute();

	// If root was a leaf, newly created nodes are also leaves
	// so prev/nex pointers must be updated
	if( root_is_leaf )
	{
		// Set newly created entries and nodes as leaves
		new_node_left->SetLeafAttribute(true);
		new_node_right->SetLeafAttribute(true);
	}
	else	// Update subtree leaves inner pointers
		UpdatePointers( new_node_left, new_node_right );

	// Update newly created nodes pointers leaves 
	dummy_first_leaf->SetNext( new_node_left->GetFirstLeaf() );
	new_node_left->GetFirstLeaf()->SetPrevious( dummy_first_leaf );
	new_node_left->GetLastLeaf()->SetNext( new_node_right->GetFirstLeaf() );
	new_node_right->GetFirstLeaf()->SetPrevious( new_node_left->GetLastLeaf() );
	new_node_right->GetLastLeaf()->SetNext( NULL );

	// Create new root and assign new entries
	CFNode* new_root = new CFNode();
	new_root->Add(new_entry_left);
	new_root->Add(new_entry_right);

	// Update new root first/last leaf
	new_root->SetFirstLeaf( new_node_left->GetFirstLeaf() );
	new_root->SetLastLeaf( new_node_right->GetLastLeaf() );
	
	// Reassign root
	root = new_root;


	// Increment nodes counter by 2 ( indeed two new nodes has been created)
	node_cnt+=2;
}


/*!
 *	\brief Split a node
 *
 *	\param node				: node to split
 *	\param closest_entry	: closest CFEntry of entry in the node entries
 *	\param entry			: (in/out) new entry to insert, will store a new entry to be inserted at parent-level in case of split
 *	\param bsplit			: (output) boolean
 *
 */
void CFTree::SplitNode( CFNode& node, CFEntry& closest_entry, CFEntry& entry, bool& bsplit )
{
	// Get a pointer to the node that is full
	CFNode* full_node = closest_entry.GetChild();
	assert( full_node != NULL );

	// Get prev/next of current full_node
	CFNode *prev = full_node->GetFirstLeaf()->GetPrevious();
	CFNode *next = full_node->GetLastLeaf()->GetNext();

	// Get the list of entries to consider (i.e. full_node entries plus the new entry)
	std::vector<CFEntry> entries;
	entries.reserve( CFTree::internal_node_max_entries + 1 );
	for( std::size_t i = 0 ; i < CFTree::internal_node_max_entries ; i++ )
		entries.push_back((full_node->GetEntries())[i]);
	entries.push_back(entry);

	// Find the farthest entry pair
	std::pair<CFEntry*, CFEntry*> farthest_pair;
	FindFarthestPair( &entries, farthest_pair );

	// Create two new nodes
	CFNode* new_node_left = new CFNode();
	CFNode* new_node_right = new CFNode();

	// Create two entries from the two new nodes
	CFEntry new_entry_left( new_node_left );
	CFEntry new_entry_right( new_node_right );

	// Redistribute old entries in the two newly created entries
	RedistributeEntries(&entries, farthest_pair, new_entry_left, new_entry_right);

	// Get leaf status of the full_node
	bool node_is_leaf = full_node->GetLeafAttribute();

	// If full_node is a leaf, newly created nodes are also leaves
	// so prev/nex pointers must be updated
	if( node_is_leaf )
	{
		// Set newly created nodes as leaves
		new_node_left->SetLeafAttribute(true);
		new_node_right->SetLeafAttribute(true);
	}
	else  // Update subtree pointers
		UpdatePointers(new_node_left, new_node_right );

	// Closest entry is divided into to new entries so it is replaced by the leftmost new entry
	std::size_t pos = node.Replace(closest_entry, new_entry_left);

	// Update new_entry_left pointers
	new_node_left->GetFirstLeaf()->SetPrevious(prev);
	if (prev != NULL)
		prev->SetNext( new_node_left->GetFirstLeaf() );

	// Check wheter the node is full
	bsplit = node.IsFull();

	// If there is room in the node, add the second entry to the node
	if( !bsplit )
	{
		// Insert the right entry after the left one
		node.Insert(new_entry_right, pos + 1);

		// Update new_entry_left and new_entry_right pointers
		new_node_left->GetLastLeaf()->SetNext( new_node_right->GetFirstLeaf() );
		new_node_right->GetFirstLeaf()->SetPrevious( new_node_left->GetLastLeaf() );
		new_node_right->GetLastLeaf()->SetNext( next );
		if (next != NULL)
			next->SetPrevious( new_node_right->GetLastLeaf() );
	}
	// else copy the second entry created into return variable 'new_entry'	
	else
	{
		// Assign new_entry_right as the new entry to considerer at parent-level
		entry = new_entry_right;

		// Update left node pointers (right node ones will be updated at parent-level)
		new_node_left->GetLastLeaf()->SetNext( next );
		if ( next != NULL)
				next->SetPrevious( new_node_left->GetLastLeaf() );
	}

	// Update node first/last leaf
	node.SetFirstLeaf( node.GetEntries()[0].GetChild()->GetFirstLeaf() );
	node.SetLastLeaf( node.GetEntries()[node.GetSize()-1].GetChild()->GetLastLeaf() );

	// Increment nodes counter by 1 (two nodes have been created, but one have "removed" with the Replace())
	node_cnt++;
}



/*!
 *	\brief Redistribute the entries into the newly created nodes
 *
 *	\param entries			: vector of entries to redistribute
 *	\param fartest_pair		: farthest pair found
 *	\param new_entry_left	: first new entry 
 *	\param new_entry_right	: second new entry
 *
 */
void CFTree::RedistributeEntries( vector<CFEntry> *entries, std::pair<CFEntry*,CFEntry*> farthest_pair, CFEntry& new_entry_left, CFEntry& new_entry_right )
{
	// Add first pair entry in new_entry_left 
	(new_entry_left.GetChild())->Add(*(farthest_pair.first));
	new_entry_left += *farthest_pair.first;
	// and the second one in the new_entry_right
	(new_entry_right.GetChild())->Add(*(farthest_pair.second));
	new_entry_right += *farthest_pair.second;

	// Process all the entries to be redistributed
	for( std::size_t i = 0 ; i < entries->size() ; i++ )
	{
		CFEntry *e = &(entries->at(i));

		// If one of the farthest pair, do not process
		if( e == farthest_pair.first || e == farthest_pair.second )
			continue;

		// Get distance to each seed
		double dist_first = CFEntry::DistD0(*e, *farthest_pair.first);
		double dist_second = CFEntry::DistD0(*e, *farthest_pair.second);

		// Assign current entry to the correct seed
		CFEntry& e_update = (dist_first < dist_second)? new_entry_left : new_entry_right;
		// Add the entry to the child node
		(e_update.GetChild())->Add(*e);
		// Update the entry
		e_update += *e;

	}
}

/*!
 *	\brief Update pointers if newly created nodes are at the next to last level
 *
 *	\param prev				: previous leaf node
 *	\param new_entry_left	: left new entry 
 *	\param new_entry_right	: right new entry
 *	\param next				: next leaf_noe
 *
 */
void CFTree::UpdatePointers( CFNode* new_node_left, CFNode* new_node_right )
{
	// Process new left entries child
	for (std::size_t i = 0 ; i < new_node_left->GetSize() - 1; i++)
	{
		new_node_left->GetEntries()[i].GetChild()->GetLastLeaf()->SetNext( new_node_left->GetEntries()[i+1].GetChild()->GetFirstLeaf() );
		new_node_left->GetEntries()[i+1].GetChild()->GetFirstLeaf()->SetPrevious( new_node_left->GetEntries()[i].GetChild()->GetLastLeaf() );
	}

	// Process new right entries child
	for (std::size_t i = 0 ; i < new_node_right->GetSize() - 1; i++)
	{
		new_node_right->GetEntries()[i].GetChild()->GetLastLeaf()->SetNext( new_node_right->GetEntries()[i+1].GetChild()->GetFirstLeaf() );
		new_node_right->GetEntries()[i+1].GetChild()->GetFirstLeaf()->SetPrevious( new_node_right->GetEntries()[i].GetChild()->GetLastLeaf() );
	}

	// Update first/last pointers newly created nodes
	new_node_left->SetFirstLeaf( new_node_left->GetEntries()[0].GetChild()->GetFirstLeaf() );
	new_node_left->SetLastLeaf( new_node_left->GetEntries()[new_node_left->GetSize()-1].GetChild()->GetLastLeaf() );
	new_node_right->SetFirstLeaf( new_node_right->GetEntries()[0].GetChild()->GetFirstLeaf() );
	new_node_right->SetLastLeaf( new_node_right->GetEntries()[new_node_right->GetSize()-1].GetChild()->GetLastLeaf() );
}

/*!
 *	\brief Rebuild tree from the existing leaf entries.
 *
 *	Rebuilding cftree is regarded as clustering, because there could be overlapped cfentries.
 *  Birch guarantees datapoints in cfentries within a range, but two data-points within a range can be separated to different cfentries
 *
 *	\param extend : if true, the size of tree reaches to memory limit, so distance threshold enlarges.
 *                  in case of both true and false, rebuilding CFTree from the existing leaves.
 *
 */
void CFTree::Rebuild( bool extend )
{
	if( extend )
	{
		// decide the next threshold
		double new_threshold = std::pow(AverageDistanceOfClosestPairLeafEntries(),2);
		threshold = threshold > new_threshold ? threshold*2 : new_threshold;
	}

	// construct a new tree by inserting all the node from the previous tree
	CFTree new_tree( dim, threshold, CFTree::internal_node_max_entries, CFTree::leaf_max_entries );
		
	CFNode* leaf = dummy_first_leaf;
	while( leaf != NULL )
	{
		for( std::size_t i = 0 ; i < leaf->GetSize() ; i++ )
			new_tree.Insert((leaf->GetEntries())[i]);

		// next leaf
		leaf = leaf->GetNext();
	}

	// Reassign tree attributes
	root = new_tree.root;
	dummy_first_leaf = new_tree.dummy_first_leaf;
	node_cnt = new_tree.node_cnt;
}


/*!
 *	\brief Return the clusters obtained
 *
 *	\return Set of clusters (one by CFEntry)
 *
 */
std::size_t CFTree::ComputeHeight()
{
	std::size_t max_height = 0;

	// Only if the root is not a leaf
	if (! root->GetLeafAttribute())
		for (std::size_t i=0; i<root->GetSize(); i++)
			max_height = std::max(max_height, (((root->GetEntries())[i]).GetChild())->GetHeight() );

	height = max_height + 1;

	return height;
}


/*!
 *	\brief Return the clusters obtained
 *
 *	\return Set of clusters (one by CFEntry)
 *
 */
std::vector<Cluster *> CFTree::GetClusters()
{
	std::size_t nleaf = 0;
	vector<Cluster *> clusters;

	// Get first real leaf
	CFNode* leaf = dummy_first_leaf->GetNext();

	// For each leaf
	while (leaf != NULL)
	{
		// For each entries within this leaf
		for (std::size_t i = 0; i < leaf->GetSize(); i++)
			clusters.push_back( ((leaf->GetEntries())[i]).GetClusters() );

		// Go to next leaf
		leaf = leaf->GetNext();
	}

	return clusters;
}


/*!
 *	\brief Return the clusters obtained without outliers
 *         An outlier is an cluster containing only one item
 *
 *	\return Set of relevant clusters
 *
 */
std::vector<Cluster *> CFTree::GetPrunedClusters()
{
	std::size_t nleaf = 0;
	vector<Cluster *> clusters;
	vector<Cluster *> tmp_clusters;


	// Get first real leaf
	CFNode* leaf = dummy_first_leaf->GetNext();

	// For each leaf
	while (leaf != NULL)
	{
		// For each entries within this leaf
		for (std::size_t i = 0; i < leaf->GetSize(); i++)
		{
			// Get clusters of the current leaf
			tmp_clusters.push_back( ((leaf->GetEntries())[i]).GetClusters() );

			// Add relevant clusters
			for (std::size_t i = 0; i < tmp_clusters.size(); i++)
				if (tmp_clusters.at(i)->GetSize() > 1 )
					clusters.push_back( tmp_clusters.at(i) );

			// Reset the tmp vector
			tmp_clusters.clear();
		}

		// Go to next leaf
		leaf = leaf->GetNext();
	}

	return clusters;
}


/*!
 *	\brief Batch mode: travel along the tree, assign initial IDs, 
 *
 */
void CFTree::AssignInitialID()
{
	// Assign root ID
	root->SetID("0");

	// Assign IDs of root entries
	for (std::size_t i=0; i<root->GetSize(); i++)
		(root->GetEntries())[i].SetID("0." + utils_to_string(i));

	// Go to next level
	if (!root->GetLeafAttribute())
		for (std::size_t i=0; i<root->GetSize(); i++)
			(root->GetEntries())[i].GetChild()->AssignInitialID();
}

/*!
 * \brief Batch mode : compute representatives for each CFEntries
 *
 */
void CFTree::ComputeRepresentatives()
{
	// Compute root entries representatives
	for (std::size_t i=0; i<root->GetSize(); i++)
		(root->GetEntries())[i].ComputeRepresentatives();
}


/*!
 *	\brief Batch mode :
 *		   (i) Travel along the tree, assign initial IDs, 
 *         (ii) create RNG for each node
 *
 */
void CFTree::CreateHierarchicalRNG()
{
	// Create RNG of root entries
	root->CreateHierarchicalRNG();

	// Go to next level
	if (!root->GetLeafAttribute())
		for (std::size_t i=0; i<root->GetSize(); i++)
			(root->GetEntries())[i].GetChild()->CreateHierarchicalRNG();
}

/*!
 *	\brief Batch mode :
 *		   (i) Travel along the tree, assign initial IDs, 
 *         (ii) create RNG for each level node entries
 *
 */
void CFTree::CreateMultilevelRNG()
{
	// Initialise list of rng and its informations. 
	rngs_nodes = std::vector<std::vector<CFEntry*>>(height + 1);
	rngs_edges = std::vector<std::vector<Edge>>(height + 1);
	rngs_nb_nodes = std::vector<std::size_t>(height + 1);
	rngs_nb_edges = std::vector<std::size_t>(height + 1);

	// Create list of nodes and entries to consider at each level
	std::list<CFNode*> current_level_nodes;
	std::list<CFNode*> next_level_nodes;
	std::list<CFEntry*> entries;
	std::size_t nb_entries;

	// Add root for the first level 
	current_level_nodes.push_back(root);

	// For each level - TODO: parallelise if possible! 
	for (std::size_t l=0; l< height + 1 ; l++)
	{
		// --------------------------------------------------------------------
		// Get list of entries of the current level. 
		// Correspond to current_level_nodes nodes entries.
		// --------------------------------------------------------------------
		// Reset number of entries
		nb_entries = 0;
		std::list<Edge>* tmp_rng; 

		// For each nodes in the current level
		for (std::size_t i=0; i<current_level_nodes.size() ;i++)
		{
			// Get current node
			CFNode* node = current_level_nodes.front();

			// Get reference of the list of entries within this node
			for (std::size_t j=0; j<node->GetSize(); j++)
				entries.push_back( &(node->GetEntries()[j]) );
				
			// Add the number of entries to considere
			nb_entries += node->GetSize();

			// Remove current node
			current_level_nodes.pop_front();
		}

		// Stop if too much objects (not relevant for visualisation)
		if (nb_entries > iRNG_MINIMUM_OBJECTS)
			break;

		// --------------------------------------------------------------------
		// Create RNG for the current level list of entries 
		// This is the code of CFNode::CreateRNG().
		// --------------------------------------------------------------------
		// Create an array of pointers to centroids
		double *pData = new double[nb_entries*dim];
		double *ptr = pData;

		// Get information of the entries of this level
		for (std::size_t i=0; i<nb_entries; i++)
		{
			ptr = pData + i*dim;
			entries.front()->GetCentroid(ptr);
			rngs_nodes[l].push_back( entries.front() );						// Push the current entry in the node list
			if (l < height + 1 )
				next_level_nodes.push_back( entries.front()->GetChild() );	// Update list of nodes of the next level 
			entries.pop_front();
		}

		// Initialize edge list structure
		tmp_rng = new list<Edge>[nb_entries];
		rngs_nb_nodes[l] = nb_entries;

		// Call RNG.cpp function
		double time = RNG::Compute_RNG(pData, nb_entries, dim, tmp_rng, &rngs_nb_edges[l]);

		// Get list of edges as a vector for export purpose
		for (std::size_t i=0; i<nb_entries; i++)
			for(list<Edge>::iterator it = tmp_rng[i].begin(); it != tmp_rng[i].end(); it++)
				rngs_edges[l].push_back(*it);

		// Reinitialize tmp_rng for next loop iteration
		tmp_rng->clear();

		// Free memory
		delete[] pData;

		// Transfer next_level nodes to current level for next loop iteration
		current_level_nodes.splice(current_level_nodes.begin(),next_level_nodes);
	}
}



/*****************************************************************************/
/*  DYNAMIC INSERTION METHODS                                                */
/*****************************************************************************/

/*!
 *	\brief Insert a CFEntry in the CFTree
 *
 *	\param entry : entry to insert in the CFTree
 *
 */
void CFTree::InsertAndUpdate(CFEntry& entry)
{
	bool bsplit;

	// Insert new entry in the root
	InsertAndUpdate(root, entry, bsplit);

	// If a split has been backpropagated up to the root, and this last one need to be split
	if( bsplit )
		SplitRootAndUpdate( entry );

	// OutOfMemory management
	std::size_t curr_mem = node_cnt * sizeof(CFNode);
	if( mem_limit > 0 && node_cnt * sizeof(CFNode) > mem_limit )
	{
		Rebuild();
	}

	return;
}


/*!
 *	\brief Insert a CFEntry in the given CFNode
 *         Since the entry is inserted in an existing tree, the node is assured to be not empty
 *
 *	\param node  : node to consider
 *	\param entry : entry to insert in the CFTree
 *
 */
void CFTree::InsertAndUpdate(CFNode* node, CFEntry& entry, bool &bsplit)
{
	// Empty node (for first insertion in the tree)
	if (node->IsEmpty())
	{
		// Add the entry to the node
		node->Add(entry);
		// Assign id
		node->ReassignID();
		// No need to split anything
		bsplit = false;
		// Return with success
		return;
	}

	// Find the closest entry among node children
	CFEntry *closest_entry = FindClosestEntry( &entry, node );

	// TODO: Here, we could add a contribution
	// If distance between the new data-point and it closest entry is too large w.r.t. a given threslhold
	// We'll rather create a new entry at this node level, instead of going deeper in the tree
	// This may avoid split operations

	// If internal node
	if (closest_entry->HasChild())
	{
		// Recursively insert the entry in the child node
		InsertAndUpdate( closest_entry->GetChild(), entry, bsplit );

		// Update CFEntry if no split
		if( !bsplit )
			*closest_entry += entry;
		// else split here
		else
			SplitNodeAndUpdate( *node, *closest_entry, entry, bsplit );

		// Update RNG
		node->UpdateRNG();
	}
	// Else leaf
	else {
		// if the leaf node can absorb the new entry
		if ( CFEntry::DistD0(*closest_entry, entry) < threshold)
		{
			// Update the ID of the inserted entry
			entry.SetID( closest_entry->GetID() + "." +  utils_to_string(closest_entry->GetObjectsNumber()) );

			// Insert entry and update the leaf in the same time
			*closest_entry += entry;

			// Update parent node RNG
			node->UpdateRNG();

			// No need to split, backpropagate this information to update at parent-level
			bsplit = false;
		}
		// else add new_entry if possible
		else if( node->GetSize() < CFTree::leaf_max_entries )
		{
			// Update the ID of the inserted entry
			entry.SetID( node->GetID() + "." +  utils_to_string(node->GetSize()) );

			// Add the entry to the node
			node->Add(entry);

			// Update RNG
			node->UpdateRNG();

			// No need to split, backpropagate this information to update at parent-level
			bsplit = false;
		}
		// if not possible, notify the split situation to the parent who will handle the actual split operation
		else	
		{
			bsplit = true;
		}
	}

	return;
}


/*!
 *	\brief Split the root
 *
 *	\param entry : new entry to insert while splitting
 *
 */
void CFTree::SplitRootAndUpdate( CFEntry& entry )
{
	// Get the list of entries to consider (i.e. root entries plus the new entry)
	std::vector<CFEntry> entries; 
	entries.reserve( CFTree::internal_node_max_entries + 1 );
	for( std::size_t i = 0 ; i < CFTree::internal_node_max_entries ; i++ )
		entries.push_back((root->GetEntries())[i]);
	entries.push_back(entry);

	// Find the farthest entry pair
	std::pair<CFEntry*, CFEntry*> farthest_pair;
	FindFarthestPair( &entries, farthest_pair );

	// Create two new nodes
	CFNode* new_node_left = new CFNode();
	CFNode* new_node_right = new CFNode();

	// Create two entries from the two new nodes
	CFEntry new_entry_left( new_node_left );
	CFEntry new_entry_right( new_node_right );

	// Redistribute old entries to new entries
	RedistributeEntries(&entries, farthest_pair, new_entry_left, new_entry_right);

	// Update representatives
	new_entry_left.UpdateRepresentatives();
	new_entry_right.UpdateRepresentatives();

	// Get leaf status of the root
	bool root_is_leaf = root->GetLeafAttribute();

	// If root was a leaf, newly created nodes are also leaves
	// so prev/nex pointers must be updated
	if( root_is_leaf )
	{
		// Set newly created entries and nodes as leaves
		new_node_left->SetLeafAttribute(true);
		new_node_right->SetLeafAttribute(true);
	}
	else	// Update subtree leaves inner pointers
		UpdatePointers( new_node_left, new_node_right );

	// Update newly created nodes pointers leaves 
	dummy_first_leaf->SetNext( new_node_left->GetFirstLeaf() );
	new_node_left->GetFirstLeaf()->SetPrevious( dummy_first_leaf );
	new_node_left->GetLastLeaf()->SetNext( new_node_right->GetFirstLeaf() );
	new_node_right->GetFirstLeaf()->SetPrevious( new_node_left->GetLastLeaf() );
	new_node_right->GetLastLeaf()->SetNext( NULL );

	// Create new root and assign new entries
	CFNode* new_root = new CFNode();
	new_root->Add(new_entry_left);
	new_root->Add(new_entry_right);

	// Update new root first/last leaf
	new_root->SetFirstLeaf( new_node_left->GetFirstLeaf() );
	new_root->SetLastLeaf( new_node_right->GetLastLeaf() );
	
	// Reassign root
	root = new_root;

	// Update ID
	root->SetID("0");
	root->ReassignID();

	// Update RNG
	root->UpdateRNG();

	// Increment nodes counter by 2 ( indeed two new nodes has been created)
	node_cnt+=2;
}


/*!
 *	\brief Split a node
 *
 *	\param node				: node to split
 *	\param closest_entry	: closest CFEntry of entry in the node entries
 *	\param entry			: (in/out) new entry to insert, will store a new entry to be inserted at parent-level in case of split
 *	\param bsplit			: (output) boolean
 *
 */
void CFTree::SplitNodeAndUpdate( CFNode& node, CFEntry& closest_entry, CFEntry& entry, bool& bsplit )
{
	// Get a pointer to the node that is full
	CFNode* full_node = closest_entry.GetChild();
	assert( full_node != NULL );

	// Get prev/next of current full_node
	CFNode *prev = full_node->GetFirstLeaf()->GetPrevious();
	CFNode *next = full_node->GetLastLeaf()->GetNext();

	// Get the list of entries to consider (i.e. full_node entries plus the new entry)
	std::vector<CFEntry> entries;
	entries.reserve( CFTree::internal_node_max_entries + 1 );
	for( std::size_t i = 0 ; i < CFTree::internal_node_max_entries ; i++ )
		entries.push_back((full_node->GetEntries())[i]);
	entries.push_back(entry);

	// Find the farthest entry pair
	std::pair<CFEntry*, CFEntry*> farthest_pair;
	FindFarthestPair( &entries, farthest_pair );

	// Create two new nodes
	CFNode* new_node_left = new CFNode();
	CFNode* new_node_right = new CFNode();

	// Create two entries from the two new nodes
	CFEntry new_entry_left( new_node_left );
	CFEntry new_entry_right( new_node_right );

	// Redistribute old entries in the two newly created entries
	RedistributeEntries(&entries, farthest_pair, new_entry_left, new_entry_right);

	// Update new nodes rng
	new_node_left->UpdateRNG();
	new_node_right->UpdateRNG();

	// Update new entries representatives
	new_entry_left.UpdateRepresentatives();
	new_entry_right.UpdateRepresentatives();

	// Get leaf status of the full_node
	bool node_is_leaf = full_node->GetLeafAttribute();

	// If full_node is a leaf, newly created nodes are also leaves
	// so prev/nex pointers must be updated
	if( node_is_leaf )
	{
		// Set newly created nodes as leaves
		new_node_left->SetLeafAttribute(true);
		new_node_right->SetLeafAttribute(true);
	}
	else  // Update subtree pointers
		UpdatePointers(new_node_left, new_node_right );

	// Closest entry is divided into to new entries so it is replaced by the leftmost new entry
	std::size_t pos = node.Replace(closest_entry, new_entry_left);

	// Update new_entry_left pointers
	new_node_left->GetFirstLeaf()->SetPrevious(prev);
	if (prev != NULL)
		prev->SetNext( new_node_left->GetFirstLeaf() );

	// Check wheter the node is full
	bsplit = node.IsFull();

	// If there is room in the node, add the second entry to the node
	if( !bsplit )
	{
		// Insert the right entry after the left one
		node.Insert(new_entry_right, pos + 1);

		// Update new_entry_left and new_entry_right pointers
		new_node_left->GetLastLeaf()->SetNext( new_node_right->GetFirstLeaf() );
		new_node_right->GetFirstLeaf()->SetPrevious( new_node_left->GetLastLeaf() );
		new_node_right->GetLastLeaf()->SetNext( next );
		if (next != NULL)
			next->SetPrevious( new_node_right->GetLastLeaf() );
	}
	// else copy the second entry created into return variable 'new_entry'	
	else
	{
		// Assign new_entry_right as the new entry to considerer at parent-level
		entry = new_entry_right;

		// Update left node pointers (right node ones will be updated at parent-level)
		new_node_left->GetLastLeaf()->SetNext( next );
		if ( next != NULL)
				next->SetPrevious( new_node_left->GetLastLeaf() );
	}

	// Update node first/last leaf
	node.SetFirstLeaf( node.GetEntries()[0].GetChild()->GetFirstLeaf() );
	node.SetLastLeaf( node.GetEntries()[node.GetSize()-1].GetChild()->GetLastLeaf() );

	// Increment nodes counter by 1 (two nodes have been created, but one have "removed" with the Replace())
	node_cnt++;

	// Update ID
	node.ReassignID();

	// Update RNG
	node.UpdateRNG();
}



/*****************************************************************************/
/*  DISTANCE RELATED METHODS                                                 */
/*****************************************************************************/

/*!
 *	\brief Find closest CFEntry of an entry among the node children
 *
 *	\param entry : entry to insert
 *	\param node	 : node to inspect
 *
 */
CFEntry* CFTree::FindClosestEntry( CFEntry* entry, CFNode* node )
{
	// Get number of entry
	int size_ = node->GetSize();

	// Get pointer to entries array
	CFEntry* begin = node->GetEntries();

	// Get closest entries
	std::size_t index = 0;
	double tmp_dist = 0.0;
	double min_dist = std::numeric_limits<double>::infinity();
	for (std::size_t i = 0; i < node->GetSize(); i++)
	{
		tmp_dist = CFEntry::DistD0(*entry, begin[i]);
		if (tmp_dist < min_dist)
		{
			min_dist = tmp_dist;
			index = i;
		}
	}

	return &((node->GetEntries())[index]);
}


/*!
 *	\brief Find fartest pair of entries 
 *
 *	\param entries		: list of entries to consider
 *	\param fartest_pair	: (output) farthest pair found
 *
 */
void CFTree::FindFarthestPair( std::vector<CFEntry> *entries, std::pair<CFEntry*, CFEntry*>& farthest_pair)
{
	assert( entries->size() >= 2 );

	double max_dist = -1.0;
	for( std::size_t i = 0 ; i < entries->size() - 1 ; i++ )
	{
		for( std::size_t j = i+1 ; j < entries->size() ; j++ )
		{
			double dist = CFEntry::DistD0( entries->at(i), entries->at(j) );
			if( dist > max_dist )
			{
				max_dist = dist;
				farthest_pair.first = &(entries->at(i));
				farthest_pair.second = &(entries->at(j));
			}
		}
	}

	return;
}


/*!
 *	\brief Get average distance of closest entries pair in leaves  
 *         for new threshold computation
 *
 */
double CFTree::AverageDistanceOfClosestPairLeafEntries()
{
	std::size_t total_n = 0;
	double	total_d = 0.0;
	double	dist;

	// Determine new threshold
	CFNode* leaf = dummy_first_leaf;
	while( leaf != NULL )
	{
		if( leaf->GetSize() >= 2 )
		{
			std::vector<double> min_dists( leaf->GetSize(), (std::numeric_limits<double>::max)() );
			for( std::size_t i = 0 ; i < leaf->GetSize() - 1 ; i++ )
			{
				for( std::size_t j = i+1 ; j < leaf->GetSize() ; j++ )
				{
					dist = CFEntry::DistD0( (leaf->GetEntries())[i], (leaf->GetEntries())[j] ); 
					dist = dist >= 0.0 ? sqrt(dist) : 0.0;
					if( min_dists[i] > dist )	min_dists[i] = dist;
					if( min_dists[j] > dist )	min_dists[j] = dist;
				}
			}
			for( std::size_t i = 0 ; i < leaf->GetSize() ; i++ )
				total_d += min_dists[i];
			total_n += leaf->GetSize();
		}

		// next leaf
		leaf = leaf->GetNext();
	}
	return total_d / total_n;
}

/*****************************************************************************/
/*  NOTIFICATION                                                             */
/*****************************************************************************/

/*!
 *	\brief Print the list of leaves
 *
 *
 */
void CFTree::PrintLeaves()
{
	std::size_t nleaf = 0;
	std::size_t nCluster = 0;

	// Get first real leaf
	CFNode* leaf = dummy_first_leaf->GetNext();

	// For each leaf
	while (leaf != NULL)
	{
		// Print
		std::cout << "LEAF " << nleaf++  << "\t #entries = " << leaf->GetSize() << "\t#objects=" << leaf->GetObjectsNumber() << std::endl;

		// Increment cluster number
		nCluster += leaf->GetSize();

		// Go to next leaf
		leaf = leaf->GetNext();
	}

	// Print
	std::cout << "\n  Number of clusters = " << nCluster << std::endl;
}

/*!
 *	\brief Print a CFTree structure in the console
 */
std::size_t CFTree::Print()
{
	std::size_t cnt = 0;

	std::cout << "\n<root id=" << root->GetID() << ">" << std::endl;

	for (std::size_t i=0; i<root->GetSize(); i++)
	{
		cnt += (root->GetEntries())[i].Print(1);
	}

	std::cout << "</root>" << std::endl;

	return cnt;
}

