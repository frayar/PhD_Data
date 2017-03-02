#include "main.h"
#include "irng.h"
#include "export.h"

// GOBALS
std::string g_nearestOutputFilename = "C:\\Users\\frédéric.rayar\\Desktop\\2k14\\_tests\\_irng_omp_mnist_d50.nearest";
std::string g_farthestOutputFilename = "C:\\Users\\frédéric.rayar\\Desktop\\2k14\\_tests\\_irng_omp_mnist_d50.farthest";
std::string g_candidatesOutputFilename = "C:\\Users\\frédéric.rayar\\Desktop\\2k14\\_tests\\_irng_omp_mnist_d50.candidates";
std::string g_debugOutputFilename = "C:\\Users\\frédéric.rayar\\Desktop\\2k14\\_tests\\_irng_omp_mnist_d50.debug";
bool g_DEBUG_OUTPUT = true;

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

	// Initialise nearest/farthest/candidates structures
	for(__int64 i = 0; i< iRow; i++)
	{
		nearest[i] = -1;
		farthest[i] = -1;
		sr_candidates[i] = -1;
	}

	// Add an edge between the two first node
	double d01 = my_distance(pfData, iCol, 0, 1);
	resultsCPU[0].insert(std::pair<__int64, double>(1,d01));
	resultsCPU[1].insert(std::pair<__int64, double>(0,d01));
	nearest[0] = farthest[0] = 1;
	nearest[1] = farthest[1] = 0;

	// DEBUG - test multithread - <<PARALLEL>>	
	std::cout << "<irng.cpp> Check multithread: ";
	#pragma omp parallel for
	for(int n=0; n<30; ++n)
		std::cout << n << " ";
	std::cout << std::endl << std::endl;

	// DEBUG
	if (g_DEBUG_OUTPUT)
	{
		ofstream debug_output(g_debugOutputFilename.c_str() , ios_base::app);
		debug_output << "Point_ID\tNN_Search\tSR\tNbCandidatesSR\tNbCheckedEdges\tStep1_time\tStep2_time\n";
		debug_output.close();
	}

	// Log insertion time
	ofstream timelog(timelogFilename , ios_base::app);

	// Incrementally insert the rest of the points
	for(__int64 i = 2; i < iRow; i++)
	{
		time1 = insert_point(i, pfData, iRow, iCol, threads, nearest, farthest, sr_candidates, resultsCPU);
		timelog << "Point " << i << " inserted in " << time1 << " seconds\n";

		// Monitor process
		if (i%50 == 0)
			std::cout << "Point " << i << " inserted in " << time1 << " seconds.\n\tProcessed time = " << ((clock() - time0)/1000.0) << " seconds." << std::endl;
	}

	// Close timelog
	timelog.close();

	// Get IRNG time
	time0 = clock() - time0;
	time0 = time0 / (double) 1000.0;

	// Export graph
	ExportGraph(resultsCPU, graphOutputFilename, iRow, exportMethod);

	// DEBUG: Export nearest
	ofstream nearest_output(g_nearestOutputFilename.c_str(), ios_base::app);
	for (__int64 i = 0; i<iRow; i++)
		nearest_output << i << "\t" << nearest[i] << std::endl;
	nearest_output << std::endl;
	nearest_output.close();

	// DEBUG: Export farthest
	ofstream farthest_output(g_farthestOutputFilename.c_str(), ios_base::app);
	for (__int64 i = 0; i<iRow; i++)
		farthest_output << i << "\t" << farthest[i] << std::endl;
	farthest_output << std::endl;
	farthest_output.close();

	// DEBUG: Export candidates
	ofstream candidates_output(g_candidatesOutputFilename.c_str(), ios_base::app);
	for (__int64 i = 0; i<iRow; i++)
		candidates_output << i << "\t" << sr_candidates[i] << std::endl;
	candidates_output << std::endl;
	candidates_output.close();

	/*
	// DEBUG: Export Distances
	double **distances;
	distances = new double*[iRow];
	for(__int64 i = 0; i < iRow; i++)
		distances[i] = new double[iRow];

	// Debug: Compute distances
	#pragma omp parallel for schedule(dynamic, 256) private(iDistance)
	for(__int64 idxI = 0; idxI < iRow; ++idxI){
		for(__int64  idxJ = idxI + 1; idxJ < iRow; ++idxJ){	
			iDistance = 0.0f;
			iDistance = my_distance(pfData, iCol, idxI, idxJ);
			distances[idxI][idxJ] = distances[idxJ][idxI] = iDistance;
		}
		// Set diagonal at 0
		distances[idxI][idxI] = 0.0f;
	}
	tmp = 0.0f;

	// Debug: Export distances
	ofstream distances_output(distancesOutputFilename , ios_base::app);
	for (__int64 i = 0; i<iRow; i++)
	{
		for (__int64 j = 0; j<iRow; j++)
			distances_output << distances[i][j] << "\t";
		distances_output << std::endl;
	}
	distances_output.close();

	// Debug: Free distances
	for(__int64 i = 0; i < iRow; i++)
		delete[] distances[i];
	delete[] distances;
	*/

	// Free memory
	delete[] nearest;
	delete[] sr_candidates;
	delete[] farthest;
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
	double nniDistance = 0.0f;
	double qiDistance = 0.0f;
	double *q_distances = new double[pointID];
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

	// Compute Search area radius for true neighbord
	double sr = (q_d_nn + my_distance(data, dimension, q_nn, farthest[q_nn])) * (1 + IRNG_EPSILON);

	// Compute list of points that may potientially have an adjacent edge to update
	std::vector<__int64> candidates;

	#pragma omp parallel for schedule(dynamic) shared(candidates) private(nniDistance)
	for(__int64 i = 0; i < pointID; i++){

		// Compute distance between i and q_nn
		nniDistance = my_distance(data, dimension, i, q_nn);

		if ( nniDistance < sr )
		{
			// <<PARALLEL>>
			#pragma omp critical(add_candidates)
			{
				candidates.push_back(i);
			}
		}

	}

	// Store computed number of candidates  
	sr_candidates[pointID] = candidates.size();
	// Get number of candidates inside SR
	__int64 nbCandidates = candidates.size();
	// Get an array from the vector to directly acces indices
	__int64 *candidates_array = &candidates[0];


		
	//////////////////////////////////////////////////
	// UPDATE SR - STEP 2
	// Add true rng edges of the newly inserted point
	//////////////////////////////////////////////////
	// TIME
	double time_update_step2 = clock();

	// Reset variables
	bool edge = true;
	__int64 currentI = -1;
	__int64 currentJ = -1;
	double distIJ,distIK,distJK;


	// Compute true neighbors of newly inserted point (i.e. between pointID=K and currentI <<PARALLEL>>
	#pragma omp parallel for schedule(dynamic) shared(resultsCPU, farthest) private(currentI,currentJ,edge,distIJ,distIK,distJK)
	for (__int64 i=0; i<nbCandidates; i++)												// WARNING: Assume point are added incrementally by ascending id. It should be rng.nbVertices!
	{
		currentI = candidates_array[i];
		distIK = q_distances[currentI];		// my_distance(currentI,pointID)
		edge = true;
		for (__int64 j=0; j<nbCandidates; j++)
		{
			currentJ = candidates_array[j];
			distIJ = my_distance(data, dimension, currentI, currentJ);
			distJK = q_distances[currentJ];		// my_distance(pointID,currentJ)
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
				resultsCPU[currentI].insert(std::pair<__int64, double>(pointID,distIK));
				resultsCPU[pointID].insert(std::pair<__int64, double>(currentI,distIK));
				// Update farthest[currI] 
				if (farthest[currentI] == -1)
					farthest[currentI] = pointID;
				else if (distIK > my_distance(data, dimension, currentI, farthest[currentI]))
					farthest[currentI] = pointID;
				// Update and farthest[pointID], 
				if (farthest[pointID] == -1)
					farthest[pointID] = currentI;
				else if (distIK > my_distance(data, dimension, pointID, farthest[pointID]))
						farthest[pointID] = currentI;
			}
		}
	}

	// TIME
	time_update_step2 = clock() - time_update_step2;



	//////////////////////////////////////////////////
	// UPDATE SR - STEP 1
	// Remove adjacent edges of candidates that does not verify RNG property because of the new point appearance
	//////////////////////////////////////////////////
	// TIME
	double time_update_step1 = clock();

	// Reset variables
	currentI = -1;
	currentJ = -1;
	__int64 currentK = -1;
	distIJ = 0.0f;
	distIK = 0.0f;
	distJK = 0.0f;
	__int64 nbUpdatedEdges = 0;

	// Update local edges strategy - <<PARALLEL>>
	#pragma omp parallel for schedule(dynamic) shared(resultsCPU, farthest, nbUpdatedEdges) private(currentI,currentJ,currentK,distIJ,distIK,distJK)
	for (__int64 i=0; i<nbCandidates; i++)												// WARNING: Assume point are added incrementally by ascending id. It should be rng.nbVertices!
	{
		currentI = candidates_array[i];

		for (__int64 j=i+1; j<nbCandidates; j++)
		{
			currentJ = candidates_array[j];

			// Check only if there is an edge between currentI and currentJ
			if ( resultsCPU[currentI].find(currentJ) != resultsCPU[currentI].end())
			{
				nbUpdatedEdges++;

				// Get distance between currentI and currentJ
				distIJ = my_distance(data, dimension, currentI, currentJ);

				for (__int64 k=0; k<nbCandidates; k++)
				{
					currentK = candidates_array[k];

					// Check the validity of the edges (currentI,currentJ) wrt. pointID
					distIK = my_distance(data, dimension, currentI, currentK);		
					distJK = my_distance(data, dimension, currentJ, currentK);
					if (distIK < distIJ && distJK < distIJ)
					{
						#pragma omp critical(remove_neighbor)
						{
							if (currentI == 91 || currentJ == 91)
								bool _stop = true;

							// Remove 2 edges in resultsCPU structure
							resultsCPU[currentI].erase(currentJ);
							resultsCPU[currentJ].erase(currentI);
							// Update farthest[currentI]
							if (farthest[currentI] == currentJ)
								farthest[currentI] = GetFarthest(resultsCPU[currentI]);
							// Update farthest[currentJ], 
							if (farthest[currentJ] == currentI)
								farthest[currentJ] = GetFarthest(resultsCPU[currentJ]);
						}
					}
				} // End of k loop
			}
		} // End of j loop
	}  // End of i loop

	// TIME
	time_update_step1 = clock() - time_update_step1;



	//////////////////////////////////////////////////
	// DEBUG
	// WARNING: increase insert_time, and consequently the total time (e.g for Iris 0.4s to 4.28s!)
	//////////////////////////////////////////////////
	if (g_DEBUG_OUTPUT)
	{
		ofstream debug_output(g_debugOutputFilename.c_str() , ios_base::app);
		debug_output << pointID << "\t" << time_nn <<  "\t" << sr << "\t" << candidates.size() << "\t" << nbUpdatedEdges  << "\t" << time_update_step1 << "\t" << time_update_step2 << "\n";
		debug_output.close();
	}

	// TIME
	time1 = clock() - time1;
	time1 = time1 / (double)1000.0;		// Divide by time.h CLOCK_PER_SECOND

	// Free memory
	delete[] q_distances;

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