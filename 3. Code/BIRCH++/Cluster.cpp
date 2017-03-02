/**
 *
 * BIRCH++ 
 * Implements BIRCH Clustering
 *
 * class Cluster
 * B-tree-like data structure consisting of summarized clusters
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-04-13
 * @version 1.0
 *
 */

#include "Utils.h"
#include "RNG.h"
#include "Cluster.h"

/*****************************************************************************/
/*  CONSTRUCTORS                                                             */
/*****************************************************************************/

/*!
 *
 *	\brief Constructor (default)
 *		   Initially, cluster have no elements
 *
 */
Cluster::Cluster()
{
	size = 0;
	elements = std::vector<Object *>();
	nb_edges = 0;
}


/*****************************************************************************/
/*  GETTER/SETTER                                                            */
/*****************************************************************************/
/*!
 *	\brief Indicate if the node is a leaf 
 *
 *	\return True id the node is a leaf, else False
 */
std::string Cluster::GetID()
{
	return id;
}

/*!
 *	\brief Define if the node is a leaf 
 *
 *	\param is_leaf_ : leaf boolean status
 *
 */
void Cluster::SetID(std::string id_)
{
	id =  id_;
	
	// Assign same ID to its elements
	for (std::size_t i = 0; i< elements.size(); i++)
		elements.at(i)->tree_id = id + "." + utils_to_string(i);
}

/*!
 *
 *	\brief Get cluster number of element
 *		
 *	\return Number of elements in this cluster
 *
 */
std::size_t Cluster::GetSize()
{
	return size;
}

/*!
 *
 *	\brief Get elements vector
 *		
 *	\return vector of objets
 *
 */
std::vector<Object *> Cluster::GetElements()
{
	return elements;
}

/*!
 *	\brief Get a pointer to the first leaf of the node subtree
 *
 *	\return Pointer to the first leaf of the node subtree
 */
std::map<__int64, double>* Cluster::GetRNG()
{
	return rng;
}

/*!
 *	\brief Get a pointer to the first leaf of the node subtree
 *
 *	\return Pointer to the first leaf of the node subtree
 */
std::size_t Cluster::GetNbEdges()
{
	return nb_edges;
}

/*****************************************************************************/
/*  OPERATOR OVERLOAD                                                        */
/*****************************************************************************/

/*!
 *	\brief Operator adding a cluster to the current Cluster
 *
 *	\param right_side_cluster : cluster to add
 *
 */
void Cluster::operator+=(Cluster* right_side_cluster)
{
	// Add right_site_cluster elements
	for( std::size_t i = 0 ; i < right_side_cluster->size ; i++ ) 
		elements.push_back(right_side_cluster->elements[i]);

	// Update cluster RNG
	UpdateRNG(right_side_cluster);

	// Increase  size of the cluster
	size += right_side_cluster->size;
}

/*****************************************************************************/
/*  METHODS                                                                  */
/*****************************************************************************/

/*!
 *
 *	\brief Add a new element reference to this cluster
 *		
 *	\param Reference to the new element to add
 *
 */
void Cluster::AddElement(Object* element)
{
	elements.push_back(element);
	size++;
}

/*!
 *
 *	\brief Add a set of element references to this cluster
 *		
 *	\param References to the new elements to add
 *
 */

void Cluster::AddElements(std::vector<Object *> elements_)
{
	for (std::size_t i = 0; i< elements_.size(); i++)
		elements.push_back(elements_.at(i));

	size += elements_.size();
}


/*****************************************************************************/
/*  RNG                                                                      */
/*****************************************************************************/

/*!
 *	\brief Batch mode :
 *         Create RNG of its entries centroid
 *
 */
void Cluster::CreateRNG(std::string entry_id)
{
	// Assign cluster ID
	id = entry_id;

	// Assign IDs of cluster objects
	for (std::size_t i=0; i<size; i++)
	{
		std::string tmpID = entry_id; tmpID += "."; tmpID += utils_to_string(i);
		elements[i]->tree_id = tmpID;
	}

	// Do not consider cluster with only one item
	if (size > 1)
	{
		// Get dimension of the data
		std::size_t dimension = elements[0]->dimension;

		// Create an array of pointers to centroids
		double *pData = new double[size*dimension];
		double *ptr = pData;

		// Copy data in a structure
		// TODO : Not good, we shall used the Objet.ptr directly!
		for (std::size_t i=0; i<size; i++)
			for (std::size_t j=0; j<dimension; j++)
				ptr[i*dimension+j] = elements[i]->ptr[j];

		// Initialize edge list structure
		rng = new std::map<__int64, double>[size];

		// Call RNG.cpp function
		double time = RNG::Compute_RNG(pData, size, dimension, rng, &nb_edges);

		// Free memory
		delete[] pData;
	}
}


/*!
 *	\brief Dynamic : Update RNG with multiple insertion
 *
 */
void Cluster::UpdateRNG(Cluster* new_cluster)
{
	// Get dimension of the data
	std::size_t dimension = elements[0]->dimension;
	
	// Get number of element to insert
	std::size_t new_size = new_cluster->GetSize();

	// Dynamically realocate a new RNG structure 
	std::map<__int64, double>* newRNG = new std::map<__int64, double>[size+new_size];

	// Allocate an array of pointers for the elements
	double *pData = new double[(size+new_size) *dimension];
	double *ptr = pData;

	// Copy current clusters data in a structure - TODO : Not good, we shall used the Objet.ptr directly!
	for (std::size_t i=0; i<size; i++)
		for (std::size_t j=0; j<dimension; j++)
			ptr[i*dimension+j] = elements[i]->ptr[j];
	// Insert all the new_cluster element one by one
	for(std::size_t i = 0; i<new_size; i++)
		for (std::size_t j=0; j<dimension; j++)
			ptr[(i+size)*dimension+j] = new_cluster->GetElements()[i]->ptr[j];

	// Special case, if only two elements in the cluster
	if (size + new_size == 2)
	{
		// Delete old rng
		if (nb_edges>1)
		{
			delete[] rng;
			nb_edges = 0;
		}

		// Initialise rs the newly created structure
		rng = newRNG;

		// Compute distance between the two elements
		double dist01 = RNG::EuclidianDist(pData, dimension, 0, 1);

		// Patch to avoid having a cluster with only two duplicate elements with distance 0.0
		//if (dist01 == 0.0) 
		//	dist01 = 0.00000000001;

		// Add an edge between the two elements
		rng[0].insert(std::pair<__int64, double>(1,dist01));
		rng[1].insert(std::pair<__int64, double>(0,dist01));

		// Set the number of edges to one
		nb_edges = 1;
	}
	// 2016-04-02 Special case to improve speed
	// If there is less than iRNG_MINIMUM_OBJECTS, use parallel O(n^3) RNG computation
	else if ( size + new_size < iRNG_MINIMUM_OBJECTS)
	{
		// Delete old rng
		if (nb_edges>1)
		{
			delete[] rng;
			nb_edges = 0;
		}

		// Reassign to newly created structure
		rng = newRNG;	

		// Recreate RNG
		RNG::Compute_RNG(pData, size + new_size, dimension, rng, &nb_edges);
	}
	// Default case: Incremental insertion of new_cluster elements one by one
	// TODO: We can make another improvment : insert smallest cluster in the larger one (less factorised code though)
	else
	{
		// Copy existing rng in the newly created structure
		std::copy(rng, rng + size, newRNG);

		// Delete old rng
		if (nb_edges>1)
		{
			delete[] rng;
			nb_edges = 0;
		}

		// Reassign to newly created structure
		rng = newRNG;

		// Insert new_cluster elements one by one
		for(std::size_t i = 0; i<new_size; i++)
			RNG::IA_RNG_Insert(pData, size+i+1, dimension, rng, &nb_edges);
	}

	// Free new_cluster rng that is not used anymore
	new_cluster->ClearRNG();

	// Free memory
	delete[] pData;

	// Exit function
	return;
}

/*!
 *	\brief Free memory by clearing RNG
 *
 */
void Cluster::ClearRNG()
{
	// Delete  rng
	if (nb_edges>1)
	{
		delete[] rng;
		nb_edges = 0;
	}
}
	
