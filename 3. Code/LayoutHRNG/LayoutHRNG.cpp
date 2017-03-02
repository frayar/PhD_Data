/**
 * 
 * Layout HRNG graphs
 *
 * Layout algorithms:
 *   - Stress Majorization ["Graph Drawing by Stress Majorization", Emden R. Gansner, Yehuda Koren, Stephen North (2003)]
 *   - Maxent ["A Maxent-Stress Model for Graph Layout", Emden R. Gansner, Florham Park, Yifan Hu, Stephen North (2012)]
 * 
 * Use of boost library
 * http://zenol.fr/blog/boost-property-tree/en.html
 *
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2016-08-17
 * @version 1.0
 *
 *
 * The programm:
 *  - loads a hrng file
 *  - layout each internal node graphs
 *  - update the positions in the hrng file
 */

// Include
#include "LayoutHRNG.h"

// Prototypes
void StressMajorizationGraphDrawingAlgorithm(Graph *g);

// Globals
//std::string json_input = "Wang_CLD_HRNG.json";
//std::string json_output = "Wang_CLD_HRNG_layout.json";

std::string json_input = "root.rng";
std::string json_output = "root_layout.rng";



/**
 * @brief Layout a graph defined by its nodes and edges and recursively layout the nodes's child graphs
 *
 * @param nodes Nodes of the graph
 * @param edges Edges of the graph
 *
 */
void Layout(ptree &nodes, ptree &edges) 
{
	// Create a boost graph instance
	Graph g;

	// Keep a (id, Vextex) map to fasten search
	std::map<std::string, Vertex> vertex_map;

	// Get collection of nodes
	for (ptree::iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		// Get information of the current node
		std:string id = it->second.get_child("id").data();
		double x = std::stod(it->second.get_child("x").data());
		double y = std::stod(it->second.get_child("y").data());

		// Add a vertex in the boost graph (with inital sm_positions)
		Vertex v = boost::add_vertex(VertexProperties(id,x,y,x,y), g);

		// Add in the vertex_map
		vertex_map[id] = v;
	}

	// Get collection of edges
	for (ptree::iterator it = edges.begin(); it != edges.end(); ++it)
	{
		// Get information of the current node
		string id = it->second.get_child("id").data();
		string source = it->second.get_child("source").data();
		string target = it->second.get_child("target").data();
		double weight = std::stod(it->second.get_child("weight").data());

		// Add a vertex in the boost graph
		std::pair<Edge, bool> e = boost::add_edge(vertex_map[source], vertex_map[target], EdgeProperties(id,weight) , g);
	}
	
	// Notification
	//std::cout << "Vertice: " << num_vertices(g) << std::endl;
    //std::cout << "Edges: " << num_edges(g) << std::endl;

	// Call Stress Majorization Drawing Algorithm
	StressMajorizationGraphDrawingAlgorithm(&g);

	// Assign computed positions to each node
	for (ptree::iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		string id = it->second.get_child("id").data();
		it->second.get_child("x").data() = boost::lexical_cast<std::string>((g[vertex_map[id]].sm_x));
		it->second.get_child("y").data() = boost::lexical_cast<std::string>((g[vertex_map[id]].sm_y));
	}

	// Go to the next level
	for (ptree::iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		// Try-catch to handle no children.nodes for leaves
		try 
		{
			// Get the child
			ptree &node_child = it->second.get_child("children");

			// Get the list of nodes
			ptree &child_nodes = node_child.get_child("nodes");

			// Get the list of edges
			ptree &child_edges = node_child.get_child("edges");

			// Parse the children
			if (child_nodes.size() > 1)
				Layout(child_nodes,child_edges);
		}
		catch (ptree_bad_path)
		{
			// Leaf
			//std::cout << "leaf!" << std::endl;
		}
	}
}

/**
 * @brief Main 
 *
 */
int main() 
{

	// Open json file
	ptree root;
    read_json(json_input, root);
	
	// Get root children nodes
	ptree &nodes =  root.get_child("nodes");

	// Get root children edges
	ptree &edges =  root.get_child("edges");

	// Parse root children nodes
	Layout(nodes,edges) ;

	// Output modified json (pretty=true by default)
	write_json(json_output, root);

	// Notification
	std::cout << "\nPress Enter to continue" << std::endl;
	std::cin.ignore();

	return 0;
}
