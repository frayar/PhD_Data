/**
 * 
 * Split an inner HRNG structure in several file
 *
 * Use of boost library
 * http://zenol.fr/blog/boost-property-tree/en.html
 *
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2016-08-23
 * @version 1.0
 *
 *
 * The programm:
 *  - loads a hrng file
 *  - Export each internal node in a separate rng file
 */

// Include
#include "SplitHRNG.h"

// Globals
std::string json_input = "C:\\Fred2\\HRNG\\run4_9349\\ImageNet_CLD_HRNG_9349.json";
std::string output_dir = "C:\\Fred2\\HRNG\\run4_9349\\ImageNet_CLD_HRNG_9349\\";


/**
 * @brief Layout a graph defined by its nodes and edges and recursively layout the nodes's child graphs
 *
 * @param nodes Nodes of the graph
 * @param edges Edges of the graph
 *
 */
void Export(std::string id, ptree &nodes, ptree &edges) 
{
	// Open file
	ofstream ofsm(output_dir + id + ".rng");

	// Write header
	ofsm << "{\n\t\"directed\": false,\n\t\"graph\": [],\n\t\"multigraph\": false,\n";

	// Process nodes
	ofsm << "\t\"nodes\": [\n";
	for (ptree::iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		// Get information of the current node
		std:string id = it->second.get_child("id").data();
		double x = std::stod(it->second.get_child("x").data());
		double y = std::stod(it->second.get_child("y").data());

		// Write information of the current node
		ofsm << "\t  {\n";
		ofsm << "\t\t\"id\": \"" + it->second.get_child("id").data() + "\",\n";
		ofsm << "\t\t\"label\": \"" + it->second.get_child("label").data() + "\",\n";
		ofsm << "\t\t\"x\": " + it->second.get_child("x").data() + ",\n";
		ofsm << "\t\t\"y\": " + it->second.get_child("y").data() + ",\n";
		ofsm << "\t\t\"size\": " + it->second.get_child("size").data() + ",\n";
		ofsm << "\t\t\"color\": \"" + it->second.get_child("color").data() + "\",\n";
		ofsm << "\t\t\"nb_images\": " + it->second.get_child("nb_images").data() + ",\n";
		ofsm << "\t\t\"representative\": \"" + it->second.get_child("representative").data() + "\",\n";
		ofsm << "\t\t\"near_representatives\": \"" + it->second.get_child("near_representatives").data() + "\",\n";
		ofsm << "\t\t\"far_representatives\": \"" + it->second.get_child("far_representatives").data() + "\",\n";
		try
		{
			ofsm << "\t\t\"first_leaf\": \"" + it->second.get_child("first_leaf").data() + "\",\n";
			ofsm << "\t\t\"last_leaf\": \"" + it->second.get_child("last_leaf").data() + "\",\n";
		}
		catch (ptree_bad_path)
		{	
			// Leaf - Do nothing
		}
		ofsm << "\t\t\"children\": \"" + it->second.get_child("id").data() + ".rng" + "\"\n";
		ofsm << "\t  },\n";
	}
	// Handle last comma issue
	ofsm.seekp(-3, std::ios_base::end);
	ofsm << "\n\t],\n";

	// Process edges
	ofsm << "\t\"edges\": [\n";
	for (ptree::iterator it = edges.begin(); it != edges.end(); ++it)
	{
		// Write information of the current node
		ofsm << "\t  {\n";
		ofsm << "\t\t\"id\": \"" + it->second.get_child("id").data() + "\",\n";
		ofsm << "\t\t\"source\": \"" + it->second.get_child("source").data() + "\",\n";
		ofsm << "\t\t\"target\": \"" + it->second.get_child("target").data() + "\",\n";
		ofsm << "\t\t\"weight\": \"" + it->second.get_child("weight").data() + "\"\n";
		ofsm << "\t  },\n";
	}
	// Handle last comma issue
	ofsm.seekp(-3, std::ios_base::end);
	ofsm << "\n\t]";

	// Write footer
	ofsm << "\n}";
	ofsm.close();



	// Go to the next level
	for (ptree::iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		// Try-catch to handle no children.nodes for leaves
		try 
		{
			// Get child id
			string id = it->second.get_child("id").data();

			// Get the child
			ptree &node_child = it->second.get_child("children");

			// Get the list of nodes
			ptree &child_nodes = node_child.get_child("nodes");

			// Get the list of edges
			ptree &child_edges = node_child.get_child("edges");

			// Export the children
			Export(id, child_nodes,child_edges);
		}
		catch (ptree_bad_path)
		{
			// Leaf
			std::cout << "leaf!" << std::endl;
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

	// Notification
	std::cout << "JSON read" << std::endl;
	/*
	// Get root children nodes
	ptree &nodes =  root.get_child("nodes");

	// Get root children edges
	ptree &edges =  root.get_child("edges");

	// Parse root children nodes
	Export("root", nodes,edges) ;
	*/
	// Notification
	std::cout << "\nPress Enter to continue" << std::endl;
	std::cin.ignore();

	return 0;
}
