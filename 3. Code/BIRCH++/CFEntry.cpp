/**
 *
 * BIRCH++ 
 * Implements BIRCH Clustering
 *
 * class CFEntry
 * Summarize a cluster of data-points
 * According to the paper, an CFEntry is given by a triple (N, LS, SS)
 * - N is the number of data-points in the cluster
 * - LS is the linear lin_sum of the N data-points (it is a vector)
 * - SS is the square lin_sum of the N data-points (it is a value)
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-04-13
 * @version 1.0
 *
 */

#include "CFEntry.h"

// Static attributes initialisation
std::size_t CFEntry::dim = 0;

// Custom distance function object
struct {
	bool operator()(std::pair<std::size_t,std::size_t> a, std::pair<std::size_t,std::size_t> b)
    {   
		return a.second > b.second;
    }   
} customGreaterSort;

struct {
	bool operator()(std::pair<std::size_t,double> a, std::pair<std::size_t,double> b)
    {   
		return a.second < b.second;
    }   
} customLesserSort;


/*****************************************************************************/
/*  CONSTRUCTORS                                                             */
/*****************************************************************************/

/*!
 *
 *	\brief Constructor (default)
 *
 */
CFEntry::CFEntry()
{
	n = 0;
	sq_sum = 0.0;

	lin_sum = new double[dim];
	std::fill(lin_sum, lin_sum + dim, 0);

	// Initialize pointer
	object = NULL;
	child = NULL;

	// Initialize one default cluster that will hold all object references until actual clustering is done
	clusters = new Cluster();

	// Initialize list of representatives
	near_representatives = vector<Object *>();
	far_representatives = vector<Object *>();
}

/*!
 *	\brief Constructor for a data-point.
 *         Initialize CFEntry with one data-point
 *
 *	\param item	: data-point to be inserted in this CFEntry

 *
 */
CFEntry::CFEntry(Object* dataPoint)
{
	// PATCH 05-04-2016 : Duplicate object to be stored in the tree (ImageNet batch process)

	// Duplicate objecy object
	object = dataPoint->clone();

	// Initialize object attributes
	double* pt = object->ptr;
	dim = object->dimension;

	// Initialize child
	child = NULL;

	n = 1;
	sq_sum = 0.0;
	for( std::size_t i = 0 ; i < dim ; i++ )
		sq_sum += pt[i] * pt[i];

	lin_sum = new double[dim];
	std::copy( pt, pt + dim, lin_sum );
	
	// Initialize one default cluster that will hold all object references until actual clustering is done
	clusters = new Cluster();
	clusters->AddElement(object);

	// Initialize list of representatives
	near_representatives = vector<Object *>();
	far_representatives = vector<Object *>();
	near_representatives.push_back(object);
	far_representatives.push_back(object);

}


/*!
 *	\brief Constructor for an existing node
 *         Initialize CFEntry with the existing node as child
 *
 *	\param node	: pointer to the existing node
 *
 */
CFEntry::CFEntry(CFNode* node)
{
	 n = 0;
	 sq_sum = 0.0;
	 lin_sum = new double[dim];
	 std::fill(lin_sum, lin_sum + dim, 0);

	 // Assign child to the existing node
	 child = node;

	 // Initialize object
	 object = NULL;

	 // Initialize cluster
	 clusters = NULL;

	 // Initialize list of representatives
	 near_representatives = vector<Object *>();
	 far_representatives = vector<Object *>();
}



/*****************************************************************************/
/*  OPERATOR OVERLOAD                                                        */
/*****************************************************************************/

/*!
 *	\brief Operator creating a new entry by merging an entry to the current CFEntry
 *
 *	\param right_side_entry : entry to add
 *
 *	\return New entry, merge of two CFEntries
 *
 */
CFEntry CFEntry::operator+( const CFEntry& right_side_entry )
{
	CFEntry newEntry;

	newEntry.n = n + right_side_entry.n;

	for( std::size_t i = 0 ; i < dim ; i++ )
		newEntry.lin_sum[i] = lin_sum[i] + right_side_entry.lin_sum[i];
	newEntry.sq_sum = sq_sum + right_side_entry.sq_sum;

	// Update cluster
	if (IsInLeaf())
	{
		// Add objects of the right_side_entry
		*newEntry.clusters += clusters;
		*newEntry.clusters += right_side_entry.clusters;
	}

	// Update entry representatives
	newEntry.UpdateRepresentatives();

	return newEntry;

}


/*!
 *	\brief Operator adding an entry to the current CFEntry
 *
 *	\param right_side_entry : entry to add
 *
 */
void CFEntry::operator+=( const CFEntry& right_side_entry )
{
	for( std::size_t i = 0 ; i < dim ; i++ ) 
	{
		double val = right_side_entry.lin_sum[i];
		lin_sum[i] += val;
		sq_sum += val*val;
	}

	n += right_side_entry.n;

	// Add objects of the right_side_entry
	if (IsInLeaf())
		*clusters += right_side_entry.clusters;

	// Update entry representatives
	UpdateRepresentatives();

}

/*!
 *	\brief Operator removing an entry from the current CFEntry
 *
 *	\param right_side_entry : entry to remove
 *
 */
void CFEntry::operator-=( const CFEntry& right_side_entry )
{
	for( std::size_t i = 0 ; i < dim ; i++ ) 
	{
		double val = right_side_entry.lin_sum[i];
		lin_sum[i] -= val;
		sq_sum -= val*val;
	}

	n -= right_side_entry.n;

	// TODO
	// Remove objects of the right_side_entry
}

/*****************************************************************************/
/*  GETTER/SETTER                                                            */
/*****************************************************************************/

/*!
 *	\brief Indicate if the node is a leaf 
 *
 *	\return True id the node is a leaf, else False
 */
std::string CFEntry::GetID()
{
	return id;
}


/*!
 *	\brief Define if the node is a leaf 
 *
 *	\param is_leaf_ : leaf boolean status
 *
 */
void CFEntry::SetID(std::string id_)
{
	id =  id_;
	
	// Next level
	if (child != NULL)
		child->SetID(id_);		// Assign same ID to its child
	else
		clusters->SetID(id);	// Assign same ID to its cluster
}


/*!
 *	\brief Get child node
 *
 *	\return Pointer to the child node 
 *
 */
CFNode* CFEntry::GetChild()
{ 
	return child;
}



/*!
 *	\brief Get object of this entry
 *
 *	\return Pointer to the child node 
 *
 */
Object* CFEntry::GetObject()
{ 
	return object;
}


/*!
 *	\brief Get number of data-points within this entry
 *
 *	\return Number of data-points within this entry
 *
 */
std::size_t CFEntry::GetObjectsNumber()
{
	return n;
}


/*!
 *	\brief Get entry data-points dimension
 *
 *	\return Entry data-points dimension
 *
 */
std::size_t CFEntry::GetDimension()
{
	return dim;
}


/*!
 *	\brief Set entry data-points dimension
 *
 *	\param Data-points dimension
 *
 */
void CFEntry::SetDimension(std::size_t dimension_)
{
	dim = dimension_;
}

/*!
 *	\brief Get clusters of this entry if it is a leaf
 *
 *	\return Clusters of the entry if it is a leaf
 *
 */
Cluster* CFEntry::GetClusters()
{
	return clusters;
}


/*!
 *	\brief Get Centroid of this entry
 *
 *
 *	\return Entry centroid
 *
 */
void CFEntry::GetCentroid(double *centroid)
{
	for (std::size_t i=0; i<dim; i++)
		centroid[i] = lin_sum[i]/n;
}

/*!
 *	\brief Get Representatives of this entry: nearest neighbors of the medoid
 *
 *
 *	\return vector of representative
 *
 */
std::vector<Object *> CFEntry::GetNearRepresentatives()
{
	return near_representatives;
}

/*!
 *	\brief Get Representatives of this entry: farthest neighbors of the medoid (CURE-like)
 *
 *
 *	\return vector of representative
 *
 */
std::vector<Object *> CFEntry::GetFarRepresentatives()
{
	return far_representatives;
}


/*****************************************************************************/
/*  OTHERS METHODS                                                           */
/*****************************************************************************/

/*!
 *	\brief Does this CFEntry has a child?
 *
 *	\return True if the entry has a child, else False  
 *
 */
bool CFEntry::HasChild()
{ 
	return child != NULL; 
}


/*!
 *
  * \brief Is the CFEntry is a leaf node ?
  *
  * \return	true if the entry is in a leaf, else false
  */
bool CFEntry::IsInLeaf()
{
	return child == NULL;
}


/*!
 * \brief Batch mode : compute representatives for this CFEntry subtree
 *
 */
void CFEntry::ComputeRepresentatives()
{
	// Compute next level entries representatives
	if ( child != NULL )
		for (std::size_t i=0; i<child->GetSize(); i++)
			child->GetEntries()[i].ComputeRepresentatives();

	// Compute this CFEntry representative
	UpdateRepresentatives();
}

/*!
 * \brief Dynamic mode : compute/update representatives for this CFEntry
 *
 */
void CFEntry::UpdateRepresentatives()
{
	// Reset representatives lists
	near_representatives.clear();
	far_representatives.clear();

	//-------------------------------------------------------------------------------------------------------
	// LEAF/CLUSTER REPRESENTATIVES
	//-------------------------------------------------------------------------------------------------------
	if (child == NULL)
	{
		if (n < 3)
		{
			for (std::size_t i=0; i<n; i++)
			{
				near_representatives.push_back(clusters->GetElements()[i]);
				far_representatives.push_back(clusters->GetElements()[i]);
			}
		}
		else
		{
			// Get centroid
			double *centroid = new double[dim];
			GetCentroid(centroid);
		
			// Get medoid
			double dist = 0.0;
			double min_dist = std::numeric_limits<double>::infinity();
			std::size_t nearest_object_index = -1;

			#pragma omp parallel for schedule(dynamic) shared(min_dist, nearest_object_index) private(dist)
			for(__int64 i = 0; i < n; i++)	// Here, 'n' is the number of objects in the cluster
			{
				// Compute distance between i and pointID
				dist = EuclidianDist( dim, centroid, clusters->GetElements()[i]->ptr );

				#pragma omp critical(update_nearest)
				{
					// Update query nearest neighbor
					if ( dist >= 0 && dist < min_dist )
					{
							nearest_object_index = i;
							min_dist = dist;
					}
				}
			}

			// Free memory 
			delete[] centroid;

			// Get maximum number of representatives
			std::size_t min_rep = (n < NB_REPRESENTATIVES)? n - 1 : NB_REPRESENTATIVES - 1 ;

			// Assign the k nearest neighbors of the medoid
			std::vector<std::pair<std::size_t,double>> distances = vector<std::pair<std::size_t,double>>(clusters->GetElements().size());
			#pragma omp parallel for schedule(dynamic) shared(distances)
				for (__int64 i = 0; i < clusters->GetElements().size(); i++)
					distances[i] = std::make_pair( i, EuclidianDist( dim, clusters->GetElements()[nearest_object_index]->ptr, clusters->GetElements()[i]->ptr ) );
			// Sort by ascending order
			std::sort(distances.begin(), distances.end(), customLesserSort );
			// Assign near representatives
			for (std::size_t i = 0; i < min_rep + 1; i++)
				near_representatives.push_back(clusters->GetElements()[distances[i].first]);

			// -----------------------------------------------
			// Farthest representatives (~ outliers)
			// -----------------------------------------------
			// Assign the medoid as the first representative
			far_representatives.push_back(clusters->GetElements()[nearest_object_index]);

			if (!CURE)
			{
				// Sort by descending order
				std::reverse(distances.begin(), distances.end());

				// Assign far representatives
				for (std::size_t i =0; i<min_rep + 1; i++)
					far_representatives.push_back(clusters->GetElements()[distances[i].first]);
			}
			else	// CURE-like farthest representatives computation
			{
				// Create vector of representatives index (to avoid duplicate representatives)
				std::vector<std::size_t> rep_index = std::vector<std::size_t>();
				rep_index.push_back(nearest_object_index);

				// Get rest of the representatives
				dist = 0.0;
				double max_dist = 0.0;
				std::size_t current_rep_index = nearest_object_index;
				std::size_t farthest_object_index = -1;
				for (std::size_t i = 0; i < min_rep; i++)
				{
					#pragma omp parallel for schedule(dynamic) shared(max_dist, farthest_object_index,rep_index) private(dist)
					for(__int64 j = 0; j < n; j++)	// Here, 'n' is the number of objects in the cluster
					{
						// Compute distance between i and pointID
						dist = EuclidianDist( dim, clusters->GetElements()[current_rep_index]->ptr, clusters->GetElements()[j]->ptr );

						#pragma omp critical(update_farthest)
						{
							// Update query nearest neighbor
							if ( dist >= 0 && dist > max_dist )
							{
								// Avoid duplicate representatives
								if ((find (rep_index.begin(), rep_index.end(), j) == rep_index.end()))
								{
									farthest_object_index = j;
									max_dist = dist;
								}
							}
						}
					}

					// PATCH - 26/09/2015  
					if (n==3 && farthest_object_index>n)
						farthest_object_index = n-1;	
					// PATCH - 30/10/2015 - NASA APOD
					if (farthest_object_index>n)
						farthest_object_index = n-1;

					// Push new representative
					far_representatives.push_back(clusters->GetElements()[farthest_object_index]);
					rep_index.push_back(farthest_object_index);

					// Reassign variables
					dist = 0.0;
					max_dist = 0.0;
					current_rep_index = farthest_object_index;
					farthest_object_index = -1;
				}
			}
		}
	}
	//-------------------------------------------------------------------------------------------------------
	// INTERNAL CFEntry REPRESENTATIVES
	//-------------------------------------------------------------------------------------------------------
	else
	{
		// Get number of entries in the next level
		std::size_t childNbEntry = child->GetSize();

		// Get list of next level entries cardinality vector<index,cardinal>
		std::vector<std::pair<std::size_t,std::size_t>> cardinal = vector<std::pair<std::size_t,std::size_t>>();
		for (std::size_t i = 0; i<childNbEntry; i++)
			cardinal.push_back(std::make_pair(i,child->GetEntries()[i].GetObjectsNumber()));
		// Sort by descending order
		std::sort(cardinal.begin(), cardinal.end(), customGreaterSort );

		// Loop on each sorted child entries
		std::size_t nb_rep_max = (n < NB_REPRESENTATIVES)? n : NB_REPRESENTATIVES;
		std::size_t nb_rep = 0;
		std::size_t entry_idx = 0;

		//while ( (nb_rep < nb_rep_max) && (entry_idx < childNbEntry) )
		while (nb_rep < nb_rep_max)
		{
			// Get current child entry number of objects
			std::size_t current_entry_size = child->GetEntries()[cardinal[entry_idx].first].GetObjectsNumber();

			// Get number of representative to get in this entry (depending on the entry size wrt the number of item in the child )
			std::size_t current_nb_rep_max = current_entry_size * nb_rep_max / n;
			if (current_nb_rep_max == 0) { current_nb_rep_max = 1; }

			// Add representative of the current child entry to the entry
			for (std::size_t j = 0; j<current_nb_rep_max && nb_rep < nb_rep_max ; j++)
			{
				near_representatives.push_back(child->GetEntries()[cardinal[entry_idx].first].GetNearRepresentatives()[j]);
				//far_representatives.push_back(child->GetEntries()[cardinal[entry_idx].first].GetFarRepresentatives()[j]);

				// Increment the counter of reprentative
				nb_rep++;
			}

			// Go to the next child entry
			entry_idx++;
			if (entry_idx >= childNbEntry)
				entry_idx -= childNbEntry;
		}

		// Farthest representatives (~ outliers)
		std::reverse(cardinal.begin(), cardinal.end());
		for (std::size_t i =0; i<childNbEntry; i++)
			far_representatives.push_back(child->GetEntries()[cardinal[i].first].GetFarRepresentatives()[0]);
	}
}

/*****************************************************************************/
/*  NOTIFICATION                                                             */
/*****************************************************************************/


/*!
 *	\brief Print a CFNode structure in the console
 */
std::size_t CFEntry::Print(std::size_t nb_tab)
{
	std::size_t cnt = 0;

	std::string tabulation = std::string(nb_tab, '\t');

	if (IsInLeaf())
	{
		std::cout << tabulation.c_str() << "<entry id=" << id.c_str() << " nb_items=" << n << "/>" << std::endl;
		cnt += n;
	}
	else
	{
		std::cout << tabulation.c_str() << "<entry id=" << id.c_str() << " >" << std::endl;
		cnt += child->Print(nb_tab+1);
		std::cout << tabulation.c_str() << "</entry>" << std::endl;
	}

	return cnt;
}




/*****************************************************************************/
/*  CLUSTER ATTRIBUTES                                                       */
/*****************************************************************************/

/*!
 *	\brief Centroid of the CFEntry
 *
 *	\param entry  : considered entry 
 *
 *	\return Entry centroid
 *
 */
Object CFEntry::Centroid( const CFEntry& entry)
{
	Object centroid;
	std::size_t dim_ = entry.dim;

	centroid.ptr = new double[dim_];
	for (std::size_t i=0; i<dim_; i++)
		centroid.ptr[i] = entry.lin_sum[i]/entry.n;

	centroid.dimension = dim_;
	centroid.label = "centroid";

	return centroid;
}


/*!
 *	\brief Radius of the CFEntry
 *         It is the average distance of data-points to the centroid
 *         TODO / TO_CHECK !
 *
 *	\param entry  : considered entry 
 *
 *	\return Entry radius
 *
 */
double CFEntry::Radius( const CFEntry& entry )
{
	if( entry.n <= 1 )
		return 0.0;

	double tmp0, tmp1 = 0.0;

	for (std::size_t i=0; i < dim ; i++)
	{
		tmp0 = entry.lin_sum[i] / entry.n;
		tmp1 += tmp0*tmp0;
	}

	double radius= entry.sq_sum/entry.n - tmp1;
		
	return std::sqrt(radius);
}


/*!
 *	\brief Diameter of the CFEntry
 *         It is the average pairwise distance within a cluster
 *         TODO / TO_CHECK !
 *
 *	\param entry  : considered entry 
 *
 *	\return Entry diameter
 *
 */
double CFEntry::Diameter( const CFEntry& entry )
{
	if( entry.n <= 1 )
		return 0.0;

	double temp = 0.0;

	for (std::size_t i = 0 ; i < dim ; i++)
		temp += entry.lin_sum[i]/entry.n * entry.lin_sum[i]/(entry.n - 1);

	double diameter = 2 * (entry.sq_sum/(entry.n - 1) - temp);

	return std::sqrt(diameter);
}


/*****************************************************************************/
/*  DISTANCES                                                                */
/*****************************************************************************/

/*!
 *	\brief Euclidean Distance of Centroid (D0)
 *
 *	\param left_side_entry  : lest-side entry
 *	\param right_side_entry : right-side entry
 *
 *	\return Euclidian distance of the two entries
 *
 */
double CFEntry::DistD0( const CFEntry& left_side_entry, const CFEntry& right_side_entry )
{
	double dist = 0.0;
	double tmp;

	for (std::size_t i = 0 ; i < dim ; i++) 
	{
		tmp =  left_side_entry.lin_sum[i]/left_side_entry.n - right_side_entry.lin_sum[i]/right_side_entry.n;
		dist += tmp*tmp;
	}

	return std::sqrt(dist);
}


/*!
 *	\brief Manhattan Distance of Centroid (D1)
 *
 *	\param left_side_entry  : lest-side entry
 *	\param right_side_entry : right-side entry
 *
 *	\return Manhanttan distance of the two entries
 *
 */
double CFEntry::DistD1( const CFEntry& left_side_entry, const CFEntry& right_side_entry )
{
	double dist = 0.0;
	double tmp;

	for (std::size_t i = 0 ; i < dim ; i++) 
	{
		tmp = std::abs(left_side_entry.lin_sum[i]/left_side_entry.n - right_side_entry.lin_sum[i]/right_side_entry.n);
		dist += tmp;
	}

	return dist;
}


/*!
 *	\brief Pairwise InterCluster Distance (D2)
 *
 *	\param left_side_entry  : lest-side entry
 *	\param right_side_entry : right-side entry
 *
 *	\return Pairwise InterCluster Distance of the two entries
 *
 */
double CFEntry::DistD2( const CFEntry& left_side_entry, const CFEntry& right_side_entry )
{
	double dot = 0.0;

	for(std::size_t i = 0 ; i < dim ; i++)
		dot += left_side_entry.lin_sum[i] * right_side_entry.lin_sum[i];

	double dist = ( right_side_entry.n*left_side_entry.sq_sum + left_side_entry.n*right_side_entry.sq_sum - 2*dot ) / (left_side_entry.n*right_side_entry.n);

	return std::sqrt(dist);
}


/*!
 *	\brief Pairwise IntraCluster Distance (D3)
 *
 *	\param left_side_entry  : lest-side entry
 *	\param right_side_entry : right-side entry
 *
 *	\return Pairwise IntraCluster Distance of the two entries
 *
 */
double CFEntry::DistD3( const CFEntry& left_side_entry, const CFEntry& right_side_entry )
{
	double tmp1, tmp2 = 0.0;
	std::size_t tmpn = left_side_entry.n+right_side_entry.n;

	for (std::size_t i = 0 ; i < dim ; i++)
	{
		tmp1 = left_side_entry.lin_sum[i] + right_side_entry.lin_sum[i];
		tmp2 += tmp1/tmpn * tmp1/(tmpn-1);
	}

	double dist = 2 * ((left_side_entry.sq_sum+right_side_entry.sq_sum)/(tmpn-1) - tmp2);

	return std::sqrt(dist);
}


/*!
 *	\brief Variance Increase Distance (D3)
 *
 *	\param left_side_entry  : lest-side entry
 *	\param right_side_entry : right-side entry
 *
 *	\return Variance Increase Distance of the two entries
 *
 */
double CFEntry::DistD4( const CFEntry& left_side_entry, const CFEntry& right_side_entry )
{
	//TODO CFEntry::D4

	return 0.0;
}

/*!
 *	\brief Euclidian distance between two numerical data-point
 *
 *	\param left_side_data_point  : lest-side data-point
 *	\param right_side_data_point : right-side data-point
 *
 *	\return Variance Increase Distance of the two entries
 *
 */
double CFEntry::EuclidianDist( const std::size_t dimension, const double* left_side_data_point, const double* right_side_data_point )
{
	double tmp = 0.0;
	double dist = 0.0;

	for(std::size_t i=0; i<dimension; i++)
	{
		tmp = left_side_data_point[i] - right_side_data_point[i];
		dist += tmp*tmp;
	}
	
	return sqrt(dist);
}
