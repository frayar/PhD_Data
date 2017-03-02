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

#include "RNG.h"


#pragma region RNG

/*****************************************************************************/
/*  RNG                                                                      */
/*****************************************************************************/

/*!
 *
 *	\brief Compute RNG with brute-force algorithm O(n^3)
 *
 * @param pfData iRow x iCol array  containing data : [ data_{0}[0..dim-1] , ... , data_{n-1}[0..dim-1] ]
 * @param nbdata Number of data considered
 * @param dimension Dimension of data
 * @param edges (out) List of edges of the computed RNG
 * @param nb_edges (out) Number of edges of the computed RNG
 *
 * @return Processing time of RNG computation algorithm in seconds 
 *
 */
double RNG::Compute_RNG(double *pData, int nbData, int dimension, list<Edge> *edges, std::size_t *nb_edges)
{

	// Attributes
	double time0;
	double **distance;
	double distIJ,distIK,distJK;
	double iDistance = 0.0f;
	double temp = 0.0f;
	bool isCanUse;
	Edge edge;
	std::size_t cpt_edges = *nb_edges;

	// Specifies the number of threads used by default in subsequent parallel sections
	omp_set_num_threads( omp_get_max_threads() );

	// Initialize similarity matrix structure
	distance = new double*[nbData];
	for(int idxI = 0; idxI < nbData; ++idxI)
		distance[idxI] = new double[nbData];

	// Set clock
	time0 = clock();

	// STEP 1: Compute similarity matrix
	#pragma omp parallel for schedule(dynamic, 256) private(iDistance, temp)
	for(int idxI = 0; idxI < nbData; idxI++){
		for(int idxJ = idxI + 1; idxJ < nbData; idxJ++){	
			iDistance = 0.0f;
			for(int counter = 0, idxIi = idxI * dimension, idxJj = idxJ * dimension; counter < dimension; counter++, idxIi++, idxJj++){
				temp = pData[idxIi] - pData[idxJj];
				temp = temp * temp;
				iDistance += temp;
			}
			distance[idxI][idxJ] = distance[idxJ][idxI] = sqrt(iDistance);
		}

		// Set diagonal at 0
		distance[idxI][idxI] = 0.0f;
	}	

	// STEP 2: for each(i,j) isCanUse if edge(i,j) exist in RNG
	#pragma omp parallel for schedule(dynamic, 256) private(edge, isCanUse,distIJ,distIK,distJK) shared(cpt_edges,edges)
	for(int idxI = 0; idxI < nbData; idxI++){		
		for(int idxJ = idxI + 1; idxJ < nbData; idxJ++){
			isCanUse = true;
			distIJ = distance[idxI][idxJ];	
			for(int idxK = 0; idxK < nbData; idxK++){
				distIK = distance[idxI][idxK];
				distJK = distance[idxJ][idxK];
				if(distIK<distIJ&&distJK<distIJ)
				{
					isCanUse = false;
					break;
				}
			}
			if(isCanUse)
			{
				// Create an edge
				edge.src = idxI;
				edge.tgt  = idxJ;
				edge.weight = distIJ;
				// Add it to the list
				edges[idxI].push_back(edge);
				// Increment counter
				cpt_edges++;
			}
		}
	}	

	// Assign new number of edges
	*nb_edges = cpt_edges;

	// Get RNG time
	time0 = clock() - time0;

	// Free memory
	for(int idxI = 0; idxI < nbData; ++idxI)
		delete[] distance[idxI];
	delete[] distance;

	return time0;
}

/*!
 *
 *	\brief BATCH: Compute RNG with brute-force algorithm O(n^3)
 *
 * @param pfData iRow x iCol array  containing data : [ data_{0}[0..dim-1] , ... , data_{n-1}[0..dim-1] ]
 * @param nbdata Number of data considered
 * @param dimension Dimension of data
 * @param edges (out) Map structure for the computed RNG edges
 * @param nb_edges (out) Number of edges of the computed RNG
 *
 * @return Processing time of RNG computation algorithm in seconds 
 *
 */
double RNG::Compute_RNG(double *pData, int nbData, int dimension, std::map<__int64, double> *rng, std::size_t *nb_edges)
{

	// Attributes
	double time0;
	double **distance;
	double distIJ,distIK,distJK;
	double iDistance = 0.0f;
	double temp = 0.0f;
	bool isCanUse;
	std::size_t cpt_edges = *nb_edges;

	// Specifies the number of threads used by default in subsequent parallel sections
	omp_set_num_threads( omp_get_max_threads() );

	// Initialize similarity matrix structure
	distance = new double*[nbData];
	for(int idxI = 0; idxI < nbData; ++idxI)
		distance[idxI] = new double[nbData];

	// Set clock
	time0 = clock();

	// STEP 1: Compute similarity matrix
	#pragma omp parallel for schedule(dynamic, 256) private(iDistance, temp)
	for(int idxI = 0; idxI < nbData; idxI++){
		for(int idxJ = idxI + 1; idxJ < nbData; idxJ++){	
			iDistance = 0.0f;
			for(int counter = 0, idxIi = idxI * dimension, idxJj = idxJ * dimension; counter < dimension; counter++, idxIi++, idxJj++){
				temp = pData[idxIi] - pData[idxJj];
				temp = temp * temp;
				iDistance += temp;
			}
			distance[idxI][idxJ] = distance[idxJ][idxI] = sqrt(iDistance);
		}

		// Set diagonal at 0
		distance[idxI][idxI] = 0.0f;
	}	

	// STEP 2: for each(i,j) isCanUse if edge(i,j) exist in RNG
	#pragma omp parallel for schedule(dynamic, 256) shared(rng,cpt_edges) private(isCanUse,distIJ,distIK,distJK)
	for(int idxI = 0; idxI < nbData; idxI++){		
		for(int idxJ = idxI + 1; idxJ < nbData; idxJ++){
			isCanUse = true;
			distIJ = distance[idxI][idxJ];	
			for(int idxK = 0; idxK < nbData; idxK++){
				distIK = distance[idxI][idxK];
				distJK = distance[idxJ][idxK];
				if(distIK<distIJ&&distJK<distIJ)
				{
					isCanUse = false;
					break;
				}
			}
			if(isCanUse)
			{
				#pragma omp critical(update_rng)
				{
					// Add tho edges in the adjacencies structure
					rng[idxI].insert(std::pair<__int64, double>(idxJ,distIJ));
					rng[idxJ].insert(std::pair<__int64, double>(idxI,distIJ));
					cpt_edges++;
				}
			}
		}
	}	

	// Assign new number of edges
	*nb_edges = cpt_edges;

	// Get RNG time
	time0 = clock() - time0;
	
	// Free memory
	for(int idxI = 0; idxI < nbData; ++idxI)
		delete[] distance[idxI];
	delete[] distance;

	return time0;
}

#pragma endregion RNG


#pragma region Incremental Approximative RNG :: INSERT

/*****************************************************************************/
/*  Incremental Approximative RNG :: INSERT NEW DATA POINT                   */
/*****************************************************************************/

/*!
 *
 *	\brief Incremental Approximative RNG algorithm
 *         Insert a new data point in an existing RNG
 *         The new data-point is in data[nbdata-1]
 *		   Note: This should be used only at cluster/leaf level to insert a new image
 *               Indeed, for internal node, having B child at most, the classic RNG will be faster
 *         Complexity: Linear wrt the number of nodes of the existing RNG
 *         Source: iRNG - REV21_V3 - February 2015 (C:\Users\rayar\Desktop\2k14\CBMI2015\Projects-Backup-19022015)
 *
 * \param pfData iRow x iCol array  containing data : [ data_{0}[0..dim-1] , ... , data_{n-1}[0..dim-1] ]
 * \param nbdata Number of data considered
 * \param dimension Dimension of data
 * \param adjacencies (in/out) Array that contains collection of each points neighbors
 * \param nb_edges (in/out) Number of edges in the RNG
 *
 * \return Processing time of IRNG insertion algorithm in seconds 
 *
 */
double RNG::IA_RNG_Insert(double *data, __int64 nbdata, __int64 dimension, std::map<__int64, double> *adjacencies, std::size_t *nb_edges)
{
	double time = clock();
	double qiDistance = 0.0f;
	double tmp = 0.0f;
	__int64 q_nn = -1;
	double q_d_nn = std::numeric_limits<double>::infinity();
	__int64 q_fn = -1;
	double q_d_fn = 0.0;
	std::size_t cpt_edges = *nb_edges;

	//------------------------------------------------------------------------
	// Nearest Neighbour Search
	//------------------------------------------------------------------------
	#pragma omp parallel for schedule(dynamic) private(qiDistance) shared(q_nn, q_d_nn)
	for(__int64 i = 0; i < nbdata-1; i++){

		// Compute distance between i and nbdata-1
		qiDistance = EuclidianDist(data, dimension, i, nbdata-1 );

		// <<PARALLEL>>
		#pragma omp critical(update_nearest)
		{
			// Update query nearest neighbor
			// ( qiDistance > 0 && qiDistance < q_d_nn )	// WARNING: qiDistance > 0 avoid duplicate as nearest neighbor
			if ( qiDistance < q_d_nn )				
			{
					q_nn = i;
					q_d_nn = qiDistance;
			}
		}
	}
	
	//------------------------------------------------------------------------
	// SR COMPUTATION
	//------------------------------------------------------------------------
	// Compute Search area radius
	double sr = ( q_d_nn + EuclidianDist(data, dimension, q_nn, GetFarthest(adjacencies[q_nn])) ) * (1 + IRNG_EPSILON);	

	// PATCH - 26/09/2015  
	if (nbdata==3 && sr==0)
		sr = (adjacencies[0].find(1)->second) * (1 + IRNG_EPSILON);	
	// PATCH - 30/10/2015 - NASA APOD
	if (sr == 0)
		sr = 2.0;

	double half_sr = sr / 2.0;

	// Compute list of points laying in SR - <<PARALLEL>>
	std::vector<__int64> candidates;
	std::vector<__int64> half_candidates;

	#pragma omp parallel for schedule(dynamic) shared(candidates, half_candidates)
	for (__int64 i=0; i<nbdata-1; i++)													// WARNING: Assume point are added incrementally by ascending id. It should be rng.nbVertices!
	{
		double tmp_dist = EuclidianDist(data, dimension, i, q_nn);
		if ( tmp_dist < sr)
		{
			// <<PARALLEL>>
			#pragma omp critical(add_candidates)
			{
				candidates.push_back(i);
				if ( tmp_dist < half_sr)
					half_candidates.push_back(i);
			}
		}
	}

	//------------------------------------------------------------------------
	// Patch - 04/08/2014
	// If at least NB_MAX points in the EPS_0-SR
	// no need to parse the EPS_1-SR
	//------------------------------------------------------------------------
	int NB_MAX = 100;
	bool half = half_candidates.size() > NB_MAX ;
	// Get an array from the vector to directly acces indices
	__int64 *candidates_array = half? &half_candidates[0]: &candidates[0];
	// Get number of candidates inside SR
	__int64 nbCandidates = half? half_candidates.size(): candidates.size();
	// Clear the non needed vector
	if (half)
		candidates.clear();
	else
		half_candidates.clear();


	//------------------------------------------------------------------------
	// STEP 1
	// Add true rng edges of the newly inserted point
	// Relatives neighbours are seeked among the SR candidates
	//------------------------------------------------------------------------
	// Reset variables
	bool edge = true;
	__int64 currentI = -1;
	__int64 currentJ = -1;
	double distIJ = 0.0f;
	double distIK = 0.0f;
	double distJK = 0.0f;
	std::vector<__int64> relative_neigbors;
	relative_neigbors.clear();

	// Compute true neighbors of newly inserted point (i.e. between pointID=K and currentI) <<PARALLEL>>
	#pragma omp parallel for schedule(dynamic) shared(adjacencies,relative_neigbors,cpt_edges) private(currentI,currentJ,edge,distIJ,distIK,distJK)
	for (__int64 i=0; i<nbCandidates; i++)
	{
		currentI = candidates_array[i];
		distIK = EuclidianDist(data, dimension, currentI, nbdata-1);
		edge = true;
		for (__int64 j=0; j<nbCandidates; j++)
		{
			currentJ = candidates_array[j];
			distIJ = EuclidianDist(data, dimension, currentI, currentJ);
			distJK = EuclidianDist(data, dimension, currentJ, nbdata-1);
			if (distIJ > 0 && distJK >0 && distIJ < distIK && distJK < distIK  )
			{
				edge = false;
				break;
			}
		}
		if (edge)
		{
			// <<PARALLEL>>
			#pragma omp critical(update_adjacencies)
			{	
				// Add tho edges in the adjacencies structure
				adjacencies[currentI].insert(std::pair<__int64, double>(nbdata-1,distIK));
				adjacencies[nbdata-1].insert(std::pair<__int64, double>(currentI,distIK));
				// Store relative neihbours of the inserted point
				relative_neigbors.push_back(currentI);	
				cpt_edges++;
			}
		}
	}

	//------------------------------------------------------------------------
	// STEP 2
	// Remove edges that does not verify RNG property because of the new point appearance
	//------------------------------------------------------------------------

	// ************************************************
	// STEP 2.1 - Retrieve list of edges to update
	// ************************************************
	// Reset variables
	currentI = -1;
	currentJ = -1;
	__int64 currentK = -1;
	__int64 currentL = -1;
	std::set<Edge*> edges;

	// Retrieve query relative neighbors as an array
	__int64 nb_q_neighbors = relative_neigbors.size();
	__int64 *q_neighbors_array = &relative_neigbors[0];
	
	// Parse all these neighbors
	#pragma omp parallel for schedule(dynamic) shared(edges) private(currentI,currentJ,currentK,currentL)
	for (__int64 i=0; i<nb_q_neighbors; i++)
	{
		// Get current neighbor id
		currentI = q_neighbors_array[i];

		// Get neighbors of currentI (2nd order vertex-neighbors of q)
		for(std::map<__int64, double>::iterator it1 = adjacencies[currentI].begin(); it1 != adjacencies[currentI].end(); it1++)
		{
			// Get currentJ id
			currentJ = it1->first;

			if(currentJ != nbdata-1)	// No need to recheck currentJ-pointID edge
			{

				// Add 2nd order edge neighbors
				#pragma omp critical(add_edge)
				{
					edges.insert(new Edge(currentI, currentJ, it1->second));
				}

				// o==3
				// Get neighbors of currentJ (2nd order vertex-neighbors of q)
				for(std::map<__int64, double>::iterator it2 = adjacencies[currentJ].begin(); it2 != adjacencies[currentJ].end(); it2++)
				{
					// Get currentK id
					currentK = it2->first;

					if(currentK != currentI)	// No need to recheck currentI-currentJ edge
					{
						// Add 3rd order edge neighbors
						#pragma omp critical(add_edge)
						{
							edges.insert(new Edge(currentJ, currentK, it2->second));
						}

						// o==4
						// Get neighbors of currentI (2nd order vertex-neighbors of q)
						for(std::map<__int64, double>::iterator it3 = adjacencies[currentK].begin(); it3 != adjacencies[currentK].end(); it3++)
						{
							// Get currentK id
							currentL = it3->first;

							if(currentL != currentJ)	// No need to recheck currentJ-currentK edge
							{
								// Add 4th order edge neighbors
								#pragma omp critical(add_edge)
								{
									edges.insert(new Edge(currentK, currentL, it3->second));
								}

								// o==5
								// ...
								// End o==5
							}
						}
						// End o==4
					}
				} //End o==3
			}
		} // End o==2
	} // End o==1


	// ************************************************************************************************
	// STEP 1.2 - Check is the edges to check are unvalidated (or not) by the newly inserted point q
	//*************************************************************************************************
	if (!edges.empty())
	{
		// Reset variables
		currentI = -1;
		currentJ = -1;
		distIJ=0.0f,distIK=0.0f,distJK=0.0f;

		// NOT PARRALLEL to avoid adjacencies stucture iterator issues
		for (std::set<Edge*>::iterator it=edges.begin(); it!=edges.end(); it++)
		{
			currentI = (*it)->src;
			currentJ = (*it)->tgt;
			distIJ = (*it)->weight;

			// Proceed only if the edge still exists
			if (adjacencies[currentI].find(currentJ) != adjacencies[currentI].end() )
			{
				// Check the validity of the edges (currentI,currentJ) wrt. pointID
				distIK = EuclidianDist(data, dimension, currentI, nbdata-1);
				distJK = EuclidianDist(data, dimension, currentJ, nbdata-1);
				if (distIK > 0 && distJK >0 && distIK < distIJ && distJK < distIJ)
				{
					// Remove 2 edges in adjacencies structure
					if ( (adjacencies[currentI].size()>1) && (adjacencies[currentJ].size()>1) )
					{
						adjacencies[currentI].erase(currentJ);
						adjacencies[currentJ].erase(currentI);
						cpt_edges--;
					}
				}
			}
		}
	}

	// Assign new number of edges
	*nb_edges = cpt_edges;

	// Debug
	//ofstream debug_output(g_debugOutputFilename.c_str() , ios_base::app);
	//debug_output << pointID << "\t" << time_nn <<  "\t" << sr << "\t" << nbCandidates << "\t" << nbRemovedEdges << "\t" << time_update_step1 << "\t" << time_update_step2 << "\n";
	//debug_output << nb_edges_to_update << "\n";
	//debug_output.close();

	// Free memory
	candidates.clear();
	half_candidates.clear();
	relative_neigbors.clear();
	
	// Time
	time = clock() - time;
	time = time / (double)1000.0;		// Divide by time.h CLOCK_PER_SECOND

	return time;
}

#pragma endregion Incremental Approximative RNG :: INSERT


#pragma region Incremental Approximative RNG :: REMOVE

/*****************************************************************************/
/*  Incremental Approximative RNG :: REMOVE EXISTING DATA POINT              */
/*****************************************************************************/

/*!
 *
 *	\brief Incremental Approximative RNG algorithm
 *         Remove an existing data point in an existing RNG
 *
 * \param data iRow x iCol array  containing data : [ data_{0}[0..dim-1] , ... , data_{n-1}[0..dim-1] ]
 * \param nbdata Number of data considered
 * \param dimension Dimension of data
 * \param data_idx Index of the data to be removed in the RNG
 * \param adjacencies (in/out) Array that contains collection of each points neighbors
 * \param nb_edges (in/out) Number of edges in the RNG
 *
 * \return Processing time of IRNG insertion algorithm in seconds 
 *
 */
double RNG::IA_RNG_Remove(double *data, __int64 nbdata, __int64 dimension, double *data_idx, std::map<__int64, double> *adjacencies, std::size_t *nb_edges)
{
	// DTODO - TODO
	return 1.0;
}

#pragma endregion Incremental Approximative RNG :: REMOVE


#pragma region Incremental Approximative RNG :: UPDATE

/*****************************************************************************/
/*  Incremental Approximative RNG :: REMOVE EXISTING DATA POINT              */
/*****************************************************************************/

/*!
 *
 *	\brief Incremental Approximative RNG algorithm
 *         Remove an existing data point in an existing RNG
 *
 * \param data iRow x iCol array  containing data : [ data_{0}[0..dim-1] , ... , data_{n-1}[0..dim-1] ]
 * \param nbdata Number of data considered
 * \param dimension Dimension of data
 * \param data_idx Index of the updated new data
 * \param adjacencies (in/out) Array that contains collection of each points neighbors
 * \param nb_edges (in/out) Number of edges in the RNG
 *
 * \return Processing time of IRNG insertion algorithm in seconds 
 *
 */
double RNG::IA_RNG_Update(double *data, __int64 nbdata, __int64 dimension, double *data_idx, std::map<__int64, double> *adjacencies, std::size_t *nb_edges)
{
	// DTODO - TODO
	//IA_RNG_Remove();
	//IA_RNG_Add();

	return 1.0;
}

#pragma endregion Incremental Approximative RNG :: UPDATE


#pragma region Utils functions

/**
  * Get farthest neighbor within the list of neighbors
  *
  * @param neighbors List of neighbors <neighbor_ID, distance>
  *
  * @returns id of the farthest neighbor
  *
  */
__int64 RNG::GetFarthest(std::map<__int64, double> neighbors)
{
	__int64 res = -1;
	double tmp_max = 0.0f;

	for (map<__int64, double>::const_iterator it = neighbors.begin(); it != neighbors.end(); it++)
	{
		if (it->second > tmp_max)
		{
			tmp_max = it->second;
			res = it->first;
		}
	}
	return res;
}


/**
  * Get nearest neighbor within the list of neighbors
  *
  * @param neighbors List of neighbors <neighbor_ID, distance>
  *
  * @returns id of the nearest neighbor
  *
  */
__int64 RNG::GetNearest(std::map<__int64, double> neighbors)
{
	__int64 res = -1;
	double tmp_min = std::numeric_limits<double>::infinity();

	for (map<__int64, double>::const_iterator it = neighbors.begin(); it != neighbors.end(); it++)
	{
		if (it->second < tmp_min)
		{
			tmp_min = it->second;
			res = it->first;
		}
	}
	return res;
}


/*!
 *	\brief Euclidian distance between between data[id1] and data[id2] (numerical data-point)

 * \param data iRow x iCol array  containing data : [ data_{0}[0..dim-1] , ... , data_{n-1}[0..dim-1] ]
 * \param dimension Dimension of data
 * \param id2 ID of the first data-point
 * \param id2 ID of the second data-point
 *
 * \return Euclidian distance between data[id1] and data[id2]
 *
 */
double RNG::EuclidianDist( double *data, std::size_t dimension, std::size_t id1, std::size_t id2)
{
	double tmp = 0.0;
	double dist = 0.0;

	for(std::size_t i=0; i<dimension; i++)
	{
		tmp = data[id1*dimension+i] - data[id2*dimension+i];
		dist += tmp*tmp;
	}
	
	return sqrt(dist);
}

#pragma endregion Utils functions