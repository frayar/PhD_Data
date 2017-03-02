/**
 * 
 * Stress majorization graph drawing algorithm
 *   ["Graph Drawing by Stress Majorization", Emden R. Gansner, Yehuda Koren, Stephen North (2003)]
 * 
 * Use of boost library
 * http://www.boost.org/doc/libs/1_61_0/libs/graph/doc/dijkstra_shortest_paths.html
 * http://programmingexamples.net/wiki/Boost/BGL/DijkstraComputePath
 *
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2016-08-18
 * @version 1.0
 *
 *
 * The programm:
 *  - loads a list of points (with an initial (x,y) placement
 *  - compute the shortest path matrix of these points
 *  - initialize a default weights matrix
 *  - layout each internal node graphs
 *  - update the positions in the hrng file
 *
 * TODO: This could be transformed in a class that will contains graph, shortestPathLengthMatrix and weightsMatrix as attributes
 * 
 */


// Include
#include "LayoutHRNG.h"

// Prototypes
void ComputeShortestPaths(Graph *g, double **shortestPathLengthMatrix);
void ComputeSMWeights(__int64 nb_nodes, double **shortestPathLengthMatrix, double **weightsMatrix);
void Go(Graph *g, double **shortestPathLengthMatrix, double **weightsMatrix);
void AtomicGo(Graph *g, double **shortestPathLengthMatrix, double **weightsMatrix);
double ComputeStress(Graph *g, double **shortestPathLengthMatrix, double **weightsMatrix);

// STRESS MAJORIZATION PARAMETERS
double EPSILON = 0.0001;
__int64 MAX_ITERATIONS = 1000;


/**
 * @brief Stress Majorization Graph Drawing Algorithm
 *
 * @param g Grapg to layout
 *
 */
void StressMajorizationGraphDrawingAlgorithm(Graph *g)
{
	// Get number of nodes
	__int64 nb_nodes = boost::num_vertices(*g);

	// Initialize shortest path matrix
	double **shortestPathLengthMatrix = new double*[nb_nodes];
	for(size_t i = 0; i < nb_nodes; i++)
	{
		shortestPathLengthMatrix[i] = new double[nb_nodes];
		for(size_t j = 0; j < nb_nodes; j++)
			shortestPathLengthMatrix[i][j] = 0.0;
	}

	// Compute shortest path as desired distances
	ComputeShortestPaths(g, shortestPathLengthMatrix);

	// Debug
	//for(size_t i = 0; i < nb_nodes; i++)
	//{
	//	for(size_t j = 0; j < nb_nodes; j++)
	//		std::cout << shortestPathLengthMatrix[i][j] << "\t";
	//	std::cout << std::endl;
	//}

	// Initialize weights matrix
	double **weightsMatrix = new double*[nb_nodes];
	for(size_t i = 0; i < nb_nodes; i++)
	{
		weightsMatrix[i] = new double[nb_nodes];
		for(size_t j = 0; j < nb_nodes; j++)
			weightsMatrix[i][j] = 0.0;
	}

	// Compute weights
	ComputeSMWeights(nb_nodes, shortestPathLengthMatrix, weightsMatrix);

	// Debug
	//for(size_t i = 0; i < nb_nodes; i++)
	//{
	//	for(size_t j = 0; j < nb_nodes; j++)
	//		std::cout << weightsMatrix[i][j] << "\t";
	//	std::cout << std::endl;
	//}

	// Start the algorithm
    Go(g, shortestPathLengthMatrix, weightsMatrix);
}



/**
 * @brief Compute Shortest path between each pair of node.
 *
 * @param g Graph to layout
 * @param shortestPathLengthMatrix (out) Matrix that contains the pairwise shortest path length
 *
 */
void ComputeShortestPaths(Graph *g, double **shortestPathLengthMatrix)
{
	// Get number of nodes
	__int64 nb_nodes = boost::num_vertices(*g);

	// Get nodes (pointers to first and last vertex)
	std::pair<boost::adjacency_list<>::vertex_iterator,boost::adjacency_list<>::vertex_iterator> vertices = boost::vertices(*g);

	// Compute Shortest paths 
	__int64 i = 0;
	for( ; vertices.first != vertices.second; ++vertices.first)
	{
		// Get node i
		Vertex src_node = *vertices.first;

		// Create things for Dijkstra
		std::vector<Vertex> predecessors(nb_nodes); // To store parents
		std::vector<double> distances(nb_nodes);    // To store distances

		// Get shortest path between node i and node j
		auto weights = boost::get(&EdgeProperties::weight,*g);
		boost::dijkstra_shortest_paths(*g, src_node, boost::distance_map(&distances[0]).weight_map(weights));

		// Assign shortest path length 
		for(size_t j = 0; j < nb_nodes; j++)
			shortestPathLengthMatrix[i][j] = shortestPathLengthMatrix[j][i] = distances[j];

		// Set diagonal at 0
		shortestPathLengthMatrix[i][i] = 0.0;

		// Increment node pointer
		i++;
	}	
};	



/**
 * @brief Compute stress majorisation algorithm weights
 *
 * @param nb_nodes Number of nodes of the considered graph
 * @param shortestPathLengthMatrix (in) Matrix that contains the pairwise shortest path length
 * @param weightsMatrix (out) Matrix that contains the pairwise stress majorization algorithm weights
 *
 * @todo Handle the alpha parameter
 */

void ComputeSMWeights(__int64 nb_nodes, double **shortestPathLengthMatrix, double **weightsMatrix)
{
	// Assign weights
	for(size_t i = 0; i < nb_nodes; i++)
	{
		for(size_t j = 0; j < nb_nodes; j++)
			weightsMatrix[i][j] = weightsMatrix[j][i] = 1; // / (shortestPathLengthMatrix[i][j] * shortestPathLengthMatrix[i][j]);
		// Set diagonal at 0
		weightsMatrix[i][i] = 0.0;
	}
}



/**
 * @brief Handle the stress majorization layout iterations
 *
 * @param g Graph to layout
 * @param shortestPathLengthMatrix (in) Matrix that contains the pairwise shortest path length
 * @param weightsMatrix (out) Matrix that contains the pairwise stress majorization algorithm weights
 *
 */

void Go(Graph *g, double **shortestPathLengthMatrix, double **weightsMatrix)
{
	// Initialise iteration counter
	__int64 cpt_iter = MAX_ITERATIONS;

	// Calculate stress for the first time
	double currentStress = ComputeStress(g, shortestPathLengthMatrix, weightsMatrix);
	double previousStress = 0.0;
	
	// Do-while loop
	do {
		// Process an iteration
        AtomicGo(g, shortestPathLengthMatrix, weightsMatrix);
		
		// Update stress value
		previousStress = currentStress;
		currentStress = ComputeStress(g, shortestPathLengthMatrix, weightsMatrix);

		//// Debug
		//if (cpt_iter == MAX_ITERATIONS)
		//{
		//	std::cout << "previousStress = " << previousStress << std::endl;
		//	std::cout << "currentStress = " << currentStress << std::endl;
		//}

		// Decrease iteration counter
		cpt_iter--;
		
	  } while ( ((previousStress-currentStress)/previousStress > EPSILON) || (cpt_iter > 0) );
}



/**
 * @brief Handle a stress majorization layout single iteration
 *
 * @param g Graph to layout
 * @param shortestPathLengthMatrix Matrix that contains the pairwise shortest path length
 * @param weightsMatrix Matrix that contains the pairwise stress majorization algorithm weights
 *
 */
void AtomicGo(Graph *g, double **shortestPathLengthMatrix, double **weightsMatrix)
{
	// Variables
	double dist = 0.0;
	double newX = 0.0;
	double newY = 0.0;
	double totalWeight = 0.0;

	// Get number of nodes
	__int64 nb_nodes = boost::num_vertices(*g);

	// Get list of sm coordinates
	auto nodes_sm_x = boost::get(&VertexProperties::sm_x,*g);
	auto nodes_sm_y = boost::get(&VertexProperties::sm_y,*g);

	// Get nodes (pointers to first and last vertex)
	std::pair<boost::adjacency_list<>::vertex_iterator,boost::adjacency_list<>::vertex_iterator> vp = boost::vertices(*g);

	// Compute Stress
	for( boost::adjacency_list<>::vertex_iterator itI = vp.first; itI != vp.second; ++itI)
	{
		// Get node i
		Vertex nodeI = *itI;

		// Reset variables
		dist = 0.0;
		newX = 0.0;
		newY = 0.0;
		totalWeight = 0.0;

		for(boost::adjacency_list<>::vertex_iterator itJ = vp.first; itJ != vp.second; ++itJ)
		{
			// Get node i
			Vertex nodeJ = *itJ;

			// Do not consider nodeI==nodeJ
			if (nodeI == nodeJ)
				continue;

			// Compute euclidian distance
			dist = (nodes_sm_x[nodeI] - nodes_sm_x[nodeJ])*(nodes_sm_x[nodeI] - nodes_sm_x[nodeJ]) + (nodes_sm_y[nodeI] - nodes_sm_y[nodeJ])*(nodes_sm_y[nodeI] - nodes_sm_y[nodeJ]);
			dist = std::sqrt( dist );

			// Compute n2 contribution to n1 on the X-axis
			newX += weightsMatrix[nodeI][nodeJ] * ( ( nodes_sm_x[nodeJ] + ( shortestPathLengthMatrix[nodeI][nodeJ] * (nodes_sm_x[nodeI] - nodes_sm_x[nodeJ]) / dist ) ) );
		  
			// Compute n2 contribution to n1 on the Y-axis
			newY += weightsMatrix[nodeI][nodeJ] * ( ( nodes_sm_y[nodeJ] + ( shortestPathLengthMatrix[nodeI][nodeJ] * (nodes_sm_y[nodeI] - nodes_sm_y[nodeJ]) / dist ) ) );
		  
			// Update total weight
			totalWeight += weightsMatrix[nodeI][nodeJ];

			//// Debug
			//std::cout << "dist = " << dist << std::endl;
			//std::cout << "nodeI = " << nodes_sm_x[nodeI] << "," << nodes_sm_y[nodeI] << std::endl;
			//std::cout << "nodeJ = " << nodes_sm_x[nodeJ] << "," << nodes_sm_y[nodeJ] << std::endl;
			//std::cout << "weightsMatrix[" << nodeI << "," << nodeJ << "] = " << "," << weightsMatrix[nodeI][nodeJ] << std::endl;
			//std::cout << "shortestPathLengthMatrix[" << nodeI << "," << nodeJ << "] = " << "," << shortestPathLengthMatrix[nodeI][nodeJ] << std::endl;
			//std::cout << "-------------------------------------------------------------------------------------------------------------" << std::endl;
		}

		// Update nodeI coordinates
		nodes_sm_x[nodeI] = newX / totalWeight;
		nodes_sm_y[nodeI] = newY / totalWeight;

		//// Debug
		//std::cout << "newX = " << newX << std::endl;
		//std::cout << "newY = " << newY << std::endl;
		//std::cout << "totalWeight = " << totalWeight << std::endl;
		//std::cout << std::endl << "-------------------------------------------------------------------------------------------------------------" << std::endl;


	}
}



/**
 * @brief Compute stress value of the system
 *
 * @param g Graph to layout
 * @param shortestPathLengthMatrix Matrix that contains the pairwise shortest path length
 * @param weightsMatrix  Matrix that contains the pairwise stress majorization algorithm weights
 *
 */
double ComputeStress(Graph *g, double **shortestPathLengthMatrix, double **weightsMatrix)
{
	// Variables
	double dist = 0.0;
	double stress = 0.0;

	// Get number of nodes
	__int64 nb_nodes = boost::num_vertices(*g);

	// Get list of sm coordinates
	auto nodes_sm_x = boost::get(&VertexProperties::sm_x,*g);
	auto nodes_sm_y = boost::get(&VertexProperties::sm_y,*g);

	// Get nodes (pointers to first and last vertex)
	std::pair<boost::adjacency_list<>::vertex_iterator,boost::adjacency_list<>::vertex_iterator> itI = boost::vertices(*g);

	// Compute Stress
	for( ; itI.first != itI.second; ++itI.first)
	{
		// Get node i
		Vertex nodeI = *itI.first;
		for(boost::adjacency_list<>::vertex_iterator itJ = itI.first ; itJ != itI.second; ++itJ)
		{
			// Get node i
			Vertex nodeJ = *itJ;

			// Compute euclidian distance
			dist = (nodes_sm_x[nodeI] - nodes_sm_x[nodeJ])*(nodes_sm_x[nodeI] - nodes_sm_x[nodeJ]) + (nodes_sm_y[nodeI] - nodes_sm_y[nodeJ])*(nodes_sm_y[nodeI] - nodes_sm_y[nodeJ]);
			dist = std::sqrt( dist );
			
			// Update stress
			stress += weightsMatrix[nodeI][nodeJ] * (dist - shortestPathLengthMatrix[nodeI][nodeJ]) * (dist - shortestPathLengthMatrix[nodeI][nodeJ]);
		}
	}
	
	return stress;
}
