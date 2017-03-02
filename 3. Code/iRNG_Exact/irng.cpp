#include "main.h"
#include "irng.h"
#include "export.h"

/**
 *
 * Compute incremental RNG from a given set of points
 * No prior graph, compute from scratch
 * No Presteps will be done here
 * WARNING: For tests purpose, here NBDATA is know ...
 *
 * @param pfData [iRow x iCol] 1D array containing data to add : [ data_{0}[0..dim-1] , ... , data_{n-1}[0..dim-1] ]
 * @param iRow	Number of data
 * @param iCol Dimension of data
 * @param threads Number of thread to use fo parallelisation
 * @param exportMethod Format name of the export
 * @param graphOutputFilename Graph output filename
 * @param edgesOutputFilename  Adjacencies output filename
 * @param distancesOutputFilename Distances matrix output filename
 * @param timelog Insertion time log
 *
 * @return Processing time of IRNG algorithm in seconds (export time is not taken into account) 
 *
 */
double Compute_iRNG(double *pfData, int iRow, int iCol, int threads, EXPORT exportMethod, const char* graphOutputFilename, const char* edgesOutputFilename, const char* distancesOutputFilename, const char* timelogFilename)
{
	// Attributes
	double time0 = 0.0;
	double time1 = 0.0;
	__int64 *nearest = new __int64[iRow];
	__int64 *farthest = new __int64[iRow];
	__int64 *sr_candidates = new __int64[iRow];
	std::map<__int64, double> *resultsCPU = NULL;

	// Set clock
	time0 = clock();

	// Specifies the number of threads used by default in subsequent parallel sections
	omp_set_num_threads(threads);

	// Initialise result structure
	resultsCPU = new std::map<__int64, double>[iRow];

	// Add an edge between the two first node
	double d01 = my_distance(pfData, iCol, 0, 1);
	resultsCPU[0].insert(std::pair<__int64, double>(1,d01));
	resultsCPU[1].insert(std::pair<__int64, double>(0,d01));

	// DEBUG - test multithread - <<PARALLEL>>	
	std::cout << "<irng.cpp> Check multithread: ";
	#pragma omp parallel for
	for(int n=0; n<30; ++n)
		std::cout << n << " ";
	std::cout << std::endl << std::endl;

	// Log insertion time
	ofstream timelog(timelogFilename , ios_base::app);

	// Incrementally insert the rest of the points
	for(__int64 i = 2; i < iRow; i++)
	{
		time1 = insert_point(i, pfData, iRow, iCol, threads, nearest, farthest, sr_candidates, resultsCPU);
		timelog << i << "\t" << time1 << "\n";

		// Monitor process
		if (i%200 == 0)
			std::cout << "Point " << i /* << " inserted in " << time1 << " seconds.\n\tProcessed time = " << ((clock() - time0)/1000.0) << " seconds."*/ << std::endl;
	}

	// Close timelog
	timelog.close();

	// Get IRNG time
	time0 = clock() - time0;
	time0 = time0 / (double) 1000.0;

	// Export graph
	ExportGraph(resultsCPU, graphOutputFilename, iRow, exportMethod);

	// Free memory
	delete[] resultsCPU;

	// Return iRNG computation
	return time0;
}





/**
 *
 * Insert a point in an existing RNG
 * No presteps
 *
 * @param pointID ID of the point to insert
 * @param pfData iRow x iCol array  containing data : [ data_{0}[0..dim-1] , ... , data_{n-1}[0..dim-1] ]
 * @param nbdata Number of data considered
 * @param dimension Dimension of data
 * @param threads Number of thread to use fo parallelisation
 * @param nearest Array with nearest adjacent neighbor of points[i] ID
 * @param farthest Array with farthest adjacent neighbor of points[i] ID
 * @param sr_candidates Array to store nb of candidates that are involved in the SR update
 * @param resultsCPU Array that contains collection of each points neighbors
 *
 * @return Processing time of IRNG insertion algorithm in seconds 
 *
 */
double insert_point(__int64 pointID, double *data, __int64 nbdata, __int64 dimension, int threads, __int64 *nearest, __int64 *farthest, __int64 *sr_candidates, std::map<__int64, double> *resultsCPU)
{
	// TIME
	double time1 = clock();
	double *q_distances = new double[pointID];											// WARNING: Assume point are added incrementally by ascending id. It should be rng.nbVertices!
	double *nn_distances = new double[pointID];										// WARNING: Assume point are added incrementally by ascending id. It should be rng.nbVertices!
	double nniDistance = 0.0f;
	double qiDistance = 0.0f;
	double tmp = 0.0f;
	__int64 q_nn = -1;
	double q_d_nn = std::numeric_limits<double>::infinity();
	__int64 q_fn = -1;
	double q_d_fn = 0.0;


	//////////////////////////////////////////////////
	// NN SEARCH - <<PARALLEL>>
	//////////////////////////////////////////////////
	// TIME
	double time_nn = clock();

	// Process new point
	#pragma omp parallel for schedule(dynamic) shared(q_distances,q_nn,q_d_nn) private(qiDistance)
	for(__int64 i = 0; i < pointID; i++){

		// Compute distance between i and pointID
		qiDistance = my_distance(data, dimension, i, pointID);
		q_distances[i] = qiDistance;

		#pragma omp critical(update_nearest)
		{
			// Update query nearest neighbor
			if ( qiDistance > 0 && qiDistance < q_d_nn )							// WARNING: qiDistance > 0 avoid duplicate as nearest neighbor (e.g. IRIS entries 143 and 102 (id: 142 - 101) --> chek!
			{
					q_nn = i;
					q_d_nn = qiDistance;
			}
		}
	}
	// Assign nearest
	nearest[pointID] = q_nn;

	// TIME
	time_nn = clock() - time_nn;
	

	//////////////////////////////////////////////////
	// CANDIDATES COMPUTATION - <<PARALLEL>>
	//////////////////////////////////////////////////
	// Compute list of points that may potientially have an adjacent edge to update
	std::vector<__int64> candidatesEdges;

	#pragma omp parallel for schedule(dynamic) shared(candidatesEdges) private(nniDistance)
	for(__int64 i = 0; i < pointID; i++){

		// Compute distance between i and q_nn
		nniDistance = my_distance(data, dimension, i, q_nn);
		nn_distances[i] = nniDistance;

		if ( (nniDistance - q_distances[i]) >= 0 )
		{
			// <<PARALLEL>>
			#pragma omp critical(add_candidates)
			{
				candidatesEdges.push_back(i);
			}
		}
	}

	// Add nn to candidates which edges may be updated
	candidatesEdges.push_back(q_nn);
	// Get number of candidates for edges updates
	__int64 nbCandidatesEdges = candidatesEdges.size();
	// Get an array from the vector to directly acces indices
	__int64 *candidates_edges_array = &candidatesEdges[0];


	//////////////////////////////////////////////////
	// UPDATE SR - STEP 1
	// Remove adjacent edges of candidates that does not verify RNG property because of the new point appearance
	//////////////////////////////////////////////////
	// TIME
	double time_update_step1 = clock();

	// Reset variables
	__int64 currentI = -1;
	__int64 currentJ = -1;
	double distIJ,distIK,distJK;
	__int64 nbUpdatedEdges = 0;

	std::vector<std::pair< std::pair<__int64, __int64>, double>> edges;

	// Get list of edges that may be impacted by the insertion of the new data - <<PARALLEL>>
	#pragma omp parallel for schedule(dynamic) shared(edges, nbUpdatedEdges,nn_distances,q_distances) private(currentI,currentJ,distIJ,distIK,distJK)
	for (__int64 i=0; i<nbCandidatesEdges; i++)												// WARNING: Assume point are added incrementally by ascending id. It should be rng.nbVertices!
	{
		currentI = candidates_edges_array[i];

		// Duplicate currentI neighbors
		std::map<__int64, double> neighborsI = resultsCPU[currentI];
		int nbNeighbors = neighborsI.size();

		// Loop on currentI neighborss
		for (map<__int64, double>::iterator ItI=neighborsI.begin(); ItI!=neighborsI.end(); ItI++)
		{
			currentJ = ItI->first;

			if ( (nn_distances[currentJ] <= nn_distances[currentI]) && (q_distances[currentJ] <= nn_distances[currentI]) )
			{
				#pragma omp critical(add_edge)
				{
					// Check the validity of the edges (currentI,currentJ) wrt. pointID
					nbUpdatedEdges++;
					edges.push_back(std::make_pair(std::make_pair(currentI,currentJ), ItI->second ));
				}
			}
		} // End of currentI neighbors
	}  // End of candidates


	// Check each edges
	__int64 nbEdges = edges.size();
	#pragma omp parallel for schedule(dynamic) shared(resultsCPU, edges) private(currentI,currentJ,distIJ,distIK,distJK)
	for (__int64 i=0; i<nbEdges; i++)			
	{
		currentI = edges[i].first.first;
		currentJ = edges[i].first.second;
		distIJ = edges[i].second;
		distIK = q_distances[currentI];			
		distJK = q_distances[currentJ];
		if (distIK < distIJ && distJK < distIJ)
		{
			#pragma omp critical(remove_neighbor)
			{
				// Remove 2 edges in resultsCPU structure
				resultsCPU[currentI].erase(currentJ);
				resultsCPU[currentJ].erase(currentI);
			}
		}
	}




	// TIME
	time_update_step1 = clock() - time_update_step1;

	
	//////////////////////////////////////////////////
	// UPDATE SR - STEP 2
	// Add true rng edges of the newly inserted point
	//////////////////////////////////////////////////
	// TIME
	double time_update_step2 = clock();

	// Reset variables
	bool edge = true;
	distIJ = 0.0f;
	distIK = 0.0f;
	distJK = 0.0f;


	// Compute true neighbors of newly inserted point (i.e. between pointID=K and currentI <<PARALLEL>>
	#pragma omp parallel for schedule(dynamic) shared(resultsCPU, farthest) private(currentI,currentJ,edge,distIJ,distIK,distJK)
	for (__int64 i=0; i<pointID; i++)												// WARNING: Assume point are added incrementally by ascending id. It should be rng.nbVertices!
	{
		distIK = q_distances[i];
		edge = true;
		for (__int64 j=0; j<pointID; j++)
		{
			distIJ = my_distance(data, dimension, i, j);
			distJK = q_distances[j];
			if ( distIJ < distIK && distJK < distIK  )
			{
				edge = false;
				break;
			}
		}
		if (edge)
		{
			// <<PARALLEL>>
			#pragma omp critical(add_neighbor)
			{
				// Add twho edges in the resultCPU structure
				resultsCPU[i].insert(std::pair<__int64, double>(pointID,distIK));
				resultsCPU[pointID].insert(std::pair<__int64, double>(i,distIK));
			}
		}
	}

	// TIME
	time_update_step2 = clock() - time_update_step2;

	// TIME
	time1 = clock() - time1;
	time1 = time1 / (double)1000.0;		// Divide by time.h CLOCK_PER_SECOND

	// Free memory
	delete[] q_distances;
	delete[] nn_distances;

	return time1;
}





/**
  * Get farthest neighbor within the list of neighbors
  *
  * @param neighbors List of neighbors <neighbor_ID, distance>
  *
  * @returns id of the farthest neighbor
  *
  */
__int64 GetFarthest(std::map<__int64, double> neighbors)
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
__int64 GetNearest(std::map<__int64, double> neighbors)
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





/**
 *
 * Compute distance between two point
 *
 * @param pfData iRow x iCol array  containing data : [ data_{0}[0..dim-1] , ... , data_{n-1}[0..dim-1] ]
 * @param dimension Dimension of data
 * @param p1_id	ID of the first point
 * @param p2_id ID of the second points
 *
 * @return Euclidian distance between p1 and p2
 *
 */
double my_distance(double *data, __int64 dimension, __int64 p1_id, __int64 p2_id)
{
	double iDistance = 0.0f;
	double temp = 0.0f;
	for(__int64 cpt = 0, ptr1 = p1_id * dimension, ptr2 = p2_id * dimension; cpt < dimension; cpt++, ptr1++, ptr2++){
		temp = data[ptr1] - data[ptr2];
		temp = temp * temp;
		iDistance += temp;
	}
	return sqrt(iDistance);
}