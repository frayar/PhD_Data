/**
 *
 * BIRCH++
 * Implements BIRCH Clustering
 * 
 * GraphExporter.cpp
 * Export graphs
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-05-30
 * @version 1.0
 *
 */

#include "Utils.h"
#include "BIRCH.h"
#include "CFTree.h"

/* PROTOTYPES */
std::string ExportNode_JSON(CFNode *node, std::size_t nb_tab);
std::string ExportCluster_JSON(Cluster *cluster, std::size_t nb_tab);


/*!
 *  \brief Export a tree as a multilevel graph (JSON)
 *
 *  \param tree  : Tree to export
 *
 */
void ExportTree_JSON(CFTree *tree, std::string filename)
{
	// Open file
	ofstream ofsm(filename);

	// Write header
	ofsm << "{\n\t\"directed\": false,\n\t\"graph\": [],\n\t\"multigraph\": false,\n";

	// Write root
	ofsm << ExportNode_JSON( tree->GetRoot(), 1 );

	// Write footer
	ofsm << "}";


	// Close file
	ofsm.close();
}


/*!
 *  \brief Export a node as a multilevel graph (JSON)
 *
 *  \param tree  : Node to export
 *
 */
std::string ExportNode_JSON(CFNode *node, std::size_t nb_tab)
{
	CFEntry tmp_entry;
	Edge tmp_edge;
	std::string node_id = node->GetID();
	std::size_t cpt = 0;
	std::string tabulation = std::string(nb_tab, '\t');
	int def_layout_div; 

	// Create resulting string
	std::string node_json;

	// Write node header
	node_json+= tabulation + "\"nodes\": [\n";

	// For each entries
	for(std::size_t i=0; i<node->GetSize(); i++)
	{
		// Get current entry
		tmp_entry = node->GetEntries()[i];

		// Compute an heuristic value for default layout
		def_layout_div = std::ceil( (double)node->GetSize() / 2 );

		// Get  current position of the node in a grid
		std::div_t position = std::div(i,def_layout_div);

		// Write the entry as a node in the RNG graph
		node_json+= tabulation + "  {\n";
		node_json+=	tabulation + "\t\"id\": \"n" + tmp_entry.GetID() + "\",\n";
		node_json+=	tabulation + "\t\"label\": \"n" + tmp_entry.GetID() + "\",\n";
		node_json+= tabulation + "\t\"x\": " + utils_to_string(position.rem) + ",\n";
		node_json+= tabulation + "\t\"y\": " + utils_to_string(position.quot) + ",\n";
		node_json+= tabulation + "\t\"size\": 1,\n";
		node_json+= tabulation + "\t\"color\": \"#3366CC\",\n";
		node_json+= tabulation + "\t\"nb_images\": " + utils_to_string(tmp_entry.GetObjectsNumber()) + ",\n";
		node_json+= tabulation + "\t\"representative\": \"" + tmp_entry.GetNearRepresentatives()[0]->imagepath + "\",\n";
		// Near Representatives attributes
		__int64 j = 0;
		node_json+= tabulation + "\t\"near_representatives\": \"";
		for (j = 0; j<tmp_entry.GetNearRepresentatives().size() - 1; j++)
			node_json+= tmp_entry.GetNearRepresentatives()[j]->imagepath + ",";
		node_json+= tmp_entry.GetNearRepresentatives()[j]->imagepath;
		node_json+= "\",\n";
		// Far Representatives attributes
		j = 0;
		node_json+= tabulation + "\t\"far_representatives\": \"";
		for (j = 0; j<tmp_entry.GetFarRepresentatives().size() - 1; j++)
			node_json+= tmp_entry.GetFarRepresentatives()[j]->imagepath + ",";
		node_json+= tmp_entry.GetFarRepresentatives()[j]->imagepath;
		node_json+= "\",\n";
		// Leaves pointers
		if (!node->GetLeafAttribute())
		{
			node_json+=	tabulation + "\t\"first_leaf\": \"n" + tmp_entry.GetChild()->GetFirstLeaf()->GetEntries()[0].GetID() + "\",\n";
			node_json+=	tabulation + "\t\"last_leaf\": \"n" + tmp_entry.GetChild()->GetLastLeaf()->GetEntries()[0].GetID() + "\",\n";
		}
		// Write children
		node_json+= tabulation + "\t\"children\": {\n";
		if (tmp_entry.HasChild())
			node_json += ExportNode_JSON(tmp_entry.GetChild(), nb_tab+2);
		else // Cluster
			node_json += ExportCluster_JSON(tmp_entry.GetClusters(), nb_tab+2);
		// Write end children
		node_json+= tabulation + "\t}\n";

		// Write end node
		node_json+= tabulation + "  },\n";
	}

	// Remove last comma
	node_json.pop_back();
	node_json.pop_back();
	node_json+= '\n';

	// Write edges
	if (node->GetNbEdges() > 0)
		{

		// Write node footer
		node_json+= tabulation + "],\n";

		// Write edge header
		node_json+= tabulation + "\"edges\": [\n";

		// Get list of edges
		std::list<Edge>* edges = node->GetRNG();

		// Write RNG edges
		for (std::size_t i=0; i<node->GetSize(); i++)
		{
			for(list<Edge>::iterator it = edges[i].begin(); it != edges[i].end(); it++)
			{
				// Get current edge
				tmp_edge = *it;

				// Write edge attributes
				node_json+= tabulation + "  {\n";
				node_json+= tabulation + "\t\"id\": \"e" + node_id + "." + utils_to_string(cpt) + "\",\n";
				node_json+= tabulation + "\t\"source\": \"n" + (node->GetEntries()[tmp_edge.src]).GetID() + "\",\n";
				node_json+= tabulation + "\t\"target\": \"n" + (node->GetEntries()[tmp_edge.tgt]).GetID() + "\",\n";
				node_json+= tabulation + "\t\"weight\": \"" + utils_to_string(tmp_edge.weight) + "\"\n";
				node_json+= tabulation + "  },\n";

				// Increment counter
				cpt++;
				}
		}

		// Remove last comma
		node_json.pop_back();
		node_json.pop_back();
		node_json+= '\n';

		// Write edge footer
		node_json+= tabulation + "]\n";
	}
	else
	{
		// Write node footer
		node_json+= tabulation + "]\n";

	}

	// Return resulting string
	return node_json;
}


/*!
 *  \brief Export a cluster as a multilevel graph (JSON)
 *
 *  \param cluster  : Cluster to export
 *
 */
std::string ExportCluster_JSON(Cluster *cluster, std::size_t nb_tab)
{
	Object *tmp_object;
	std::string cluster_id = cluster->GetID();
	std::size_t cpt = 0;
	std::string tabulation = std::string(nb_tab, '\t');
	int def_layout_div; 

	// Create resulting string
	std::string cluster_json;

	// Write node header
	cluster_json+= tabulation + "\"nodes\": [\n";

	// For each entries
	for(std::size_t i=0; i<cluster->GetSize(); i++)
	{
		// Get current entry
		tmp_object = cluster->GetElements()[i];

		// Compute an heuristic value for default layout
		def_layout_div = std::ceil( (double)cluster->GetSize() / 2 );

		// Get  current position of the node in a grid
		std::div_t position = std::div(i,def_layout_div);

		// Write the entry as a node in the RNG graph
		cluster_json+= tabulation + "  {\n";
		cluster_json+= tabulation + "\t\"id\": \"n" + tmp_object->tree_id + "\",\n";
		cluster_json+= tabulation + "\t\"label\": \"n" + tmp_object->tree_id + "\",\n";
		cluster_json+= tabulation + "\t\"x\": " + utils_to_string(position.rem) + ",\n";
		cluster_json+= tabulation + "\t\"y\": " + utils_to_string(position.quot) + ",\n";
		cluster_json+= tabulation + "\t\"size\": 1,\n";
		cluster_json+= tabulation + "\t\"representative\": \"" + tmp_object->imagepath +"\"\n";
		cluster_json+= tabulation + "  },\n";
	}

	// Remove last comma
	cluster_json.pop_back();
	cluster_json.pop_back();
	cluster_json+= '\n';

	// Write edges
	if (cluster->GetNbEdges() > 0)
	{
		// Write node footer
		cluster_json+= tabulation + "],\n";

		// Write edge header
		cluster_json+= tabulation + "\"edges\": [\n";

		// Get list of edges
		std::map<__int64, double>* edges = cluster->GetRNG();

		// Write RNG edges
		for (std::size_t i=0; i<cluster->GetSize(); i++)
		{
			for(std::map<__int64, double>::iterator it = edges[i].begin(); it != edges[i].end(); it++)
			{
				// Only export edge if source.id < target.id
				if (it->first > i)
				{
					// Write edge attributes
					cluster_json+= tabulation + "  {\n";
					cluster_json+= tabulation + "\t\"id\": \"e" + cluster_id + "." + utils_to_string(cpt) + "\",\n";
					cluster_json+= tabulation + "\t\"source\": \"n" + (cluster->GetElements()[i])->tree_id + "\",\n";
					cluster_json+= tabulation + "\t\"target\": \"n" + (cluster->GetElements()[it->first])->tree_id + "\",\n";
					cluster_json+= tabulation + "\t\"weight\": \"" + utils_to_string(it->second) + "\"\n";
					cluster_json+= tabulation + "  },\n";

					// Increment counter
					cpt++;
				}
			}
		}

		// Remove last comma
		cluster_json.pop_back();
		cluster_json.pop_back();
		cluster_json+= '\n';

		// Write edge footer
		cluster_json+= tabulation + "]\n";
	}
	else
	{
		// Write node footer
		cluster_json+= tabulation + "]\n";
	}

	// Return resulting string
	return cluster_json;

}


/*!
 *  \brief Export a node as a multilevel graph (JSON)
 *
 *  \param tree  : Node to export
 *
 */
void ExportLeaves_JSON(CFTree *tree, std::string filename)
{
	CFEntry tmp_entry;
	Edge tmp_edge;
	std::size_t cpt = 0;
	std::div_t position;
	std::string str_debug = "\n\n";

	// Open file
	ofstream ofsm(filename);

	// Write header
	ofsm << "{\n\t\"directed\": false,\n\t\"graph\": [],\n\t\"multigraph\": false,\n";

	// Get first real leaf
	CFNode* leaf_node = tree->GetDummyFirstLeaf()->GetNext();

	// Write node header
	ofsm << "\t\"nodes\": [\n";

	// loop on each leaf
	while (leaf_node != NULL)
	{
		// For each entries
		for(std::size_t i=0; i<leaf_node->GetSize(); i++)
		{
			// Get current entry
			tmp_entry = leaf_node->GetEntries()[i];
			
			// Debug 
			str_debug += tmp_entry.GetID() + "\n";

			// Increment counter
			cpt++;

			// Get  current position of the node in a grid
			position = std::div(cpt,50);

			// Write the entry as a node in the RNG graph
			ofsm << "\t  {\n";
			ofsm << "\t\t\"id\": \"n" + tmp_entry.GetID() + "\",\n";
			ofsm << "\t\t\"label\": \"n" + tmp_entry.GetID() + "\",\n";
			ofsm << "\t\t\"x\": " + utils_to_string(position.rem) + ",\n";	// or std::rand()%100
			ofsm << "\t\t\"y\": " + utils_to_string(position.quot) + ",\n";
			ofsm << "\t\t\"size\": 1,\n";
			ofsm << "\t\t\"color\": \"#3366CC\",\n";
			ofsm << "\t\t\"nb_images\": " + utils_to_string(tmp_entry.GetObjectsNumber()) + ",\n";
			ofsm << "\t\t\"representative\": \"" + tmp_entry.GetNearRepresentatives()[0]->imagepath + "\",\n";
			// Near Representatives attributes
			__int64 j = 0;
			ofsm << "\t\t\"near_representatives\": \"";
			for (j = 0; j<tmp_entry.GetNearRepresentatives().size() - 1; j++)
				ofsm << tmp_entry.GetNearRepresentatives()[j]->imagepath + ",";
			ofsm << tmp_entry.GetNearRepresentatives()[j]->imagepath;
			ofsm << "\",\n";
			// Far Representatives attributes
			j = 0;
			ofsm << "\t\t\"far_representatives\": \"";
			for (j = 0; j<tmp_entry.GetFarRepresentatives().size() - 1; j++)
				ofsm << tmp_entry.GetFarRepresentatives()[j]->imagepath + ",";
			ofsm << tmp_entry.GetFarRepresentatives()[j]->imagepath;
			ofsm << "\",\n";

			// Write children - Here correspond to a  Cluster
			ofsm << "\t\t\"children\": {\n";
			ofsm << ExportCluster_JSON(tmp_entry.GetClusters(), 3);
			ofsm << "\t\t}\n";

			ofsm << "\t  },\n";
		}

		/*
		// Write edges
		if (leaf_node->GetNbEdges() > 0)
		{
			// Write node footer
			ofsm <<  + "\t],\n";

			// Write edge header
			ofsm <<  + "\t\"edges\": [\n";

			// Get list of edges
			std::list<Edge>* edges = leaf_node->GetRNG();

			// Write RNG edges
			for (std::size_t i=0; i<leaf_node->GetSize(); i++)
			{
				for(list<Edge>::iterator it = edges[i].begin(); it != edges[i].end(); it++)
				{
					// Get current edge
					tmp_edge = *it;

					// Write edge attributes
					ofsm << "\t  {\n";
					ofsm << "\t\t\"id\": \"e" + leaf_node->GetID() + "." + utils_to_string(cpt) + "\",\n";
					ofsm << "\t\t\"source\": \"n" + (leaf_node->GetEntries()[tmp_edge.src]).GetID() + "\",\n";
					ofsm << "\t\t\"target\": \"n" + (leaf_node->GetEntries()[tmp_edge.tgt]).GetID() + "\",\n";
					ofsm << "\t\t\"weight\": \"" + utils_to_string(tmp_edge.weight) + "\"\n";
					ofsm << "\t  },\n";

					// Increment counter
					cpt++;
					}
			}

			// Write edge footer
			ofsm << "\t]\n";
		}
		else
		{
			// Write node footer
			ofsm << "]\n";
		}
		*/

		// Go to next leaf
		leaf_node = leaf_node->GetNext();

	} // End of while loop

	// Remove last node comma
	long pos = ofsm.tellp();
	ofsm.seekp(pos - 3);

	// Write node footer
	ofsm << "\n\t]\n";

	// Write footer
	ofsm << "}";

	// Debug
	ofsm << str_debug;

	// Close file
	ofsm.close();
}


/*!
 *  \brief Export a graph (JSON)
 *         Used for MLRNG
 *
 *  \param level	: Level of the graph (if none, empty string)
 *  \param nodes	: Nodes of the graph
 *  \param edges	: Edges of the graph
 *  \param filename : Output filename
 *
 */
void ExportGraph_JSON(std::string level, std::vector<CFEntry*> nodes, std::vector<Edge> edges, std::string filename)
{
	CFEntry *tmp_node;
	Edge tmp_edge;
	std::size_t cpt = 0;
	int def_layout_div; 
	

	// Open file
	std::ofstream ofsm(filename);

	// Write header
	ofsm << "{\n\t\"directed\": false,\n\t\"graph\": [],\n\t\"multigraph\": false,\n";

	// Write nodes
	ofsm << "\t\"nodes\": [\n";

	// For each nodes
	for(std::size_t i=0; i<nodes.size(); i++)
	{
		// Get current entry
		tmp_node = nodes.at(i);

		// Compute an heuristic value for default layout
		def_layout_div = std::ceil( (double)nodes.size() / 2 );

		// Get  current position of the node in a grid
		std::div_t position = std::div(i,def_layout_div);

		// Write the entry as a node in the RNG graph
		ofsm << "\t  {\n";
		ofsm << "\t\t\"id\": \"n" + tmp_node->GetID() + "\",\n";
		ofsm << "\t\t\"label\": \"n" + tmp_node->GetID() + "\",\n";
		ofsm << "\t\t\"x\": " + utils_to_string(position.rem) + ",\n";
		ofsm << "\t\t\"y\": " + utils_to_string(position.quot) + ",\n";
		ofsm << "\t\t\"size\": 1,\n";
		ofsm << "\t\t\"color\": \"#3366CC\",\n";
		ofsm << "\t\t\"nb_images\": " + utils_to_string(tmp_node->GetObjectsNumber()) + ",\n";
		ofsm << "\t\t\"representative\": \"" + tmp_node->GetNearRepresentatives()[0]->imagepath + "\",\n";
		// Near Representatives attributes
		__int64 j = 0;
		ofsm << "\t\t\"near_representatives\": \"";
		for (j = 0; j<tmp_node->GetNearRepresentatives().size() - 1; j++)
			ofsm << tmp_node->GetNearRepresentatives()[j]->imagepath + ",";
		ofsm << tmp_node->GetNearRepresentatives()[j]->imagepath;
		ofsm <<  "\",\n";
		// Far Representatives attributes
		j = 0;
		ofsm << "\t\t\"far_representatives\": \"";
		for (j = 0; j<tmp_node->GetFarRepresentatives().size() - 1; j++)
			ofsm <<  tmp_node->GetFarRepresentatives()[j]->imagepath + ",";
		ofsm <<  tmp_node->GetFarRepresentatives()[j]->imagepath;
		ofsm <<  "\",\n";
		// Leaves pointers
		if (!tmp_node->IsInLeaf())
		{
			ofsm << "\t\t\"first_leaf\": \"n" + tmp_node->GetChild()->GetFirstLeaf()->GetEntries()[0].GetID() + "\",\n";
			ofsm << "\t\t\"last_leaf\": \"n" + tmp_node->GetChild()->GetLastLeaf()->GetEntries()[0].GetID() + "\",\n";
		}
		// Write children to allow visualisation
		ofsm << "\t\t\"children\": {\n\t\t}\n";
		// Write end node
		ofsm << "\t  },\n";
	}

	// Remove last comma
	long pos = ofsm.tellp();
	ofsm.seekp(pos-3);
	ofsm << '\n';

	// Write edges
	if (edges.size() > 0)
	{
		// Write node footer
		ofsm << "\t],\n";

		// Write edge header
		ofsm << "\t\"edges\": [\n";

		// Write RNG edges
		for(std::size_t i=0; i<edges.size(); i++)
		{
			// Get current edge
			tmp_edge = edges.at(i);

			// Write edge attributes
			ofsm << "\t  {\n";
			ofsm << "\t\t\"id\": \"e" + level + "." + utils_to_string(cpt) + "\",\n";
			ofsm << "\t\t\"source\": \"n" + nodes.at(tmp_edge.src)->GetID() + "\",\n";
			ofsm << "\t\t\"target\": \"n" + nodes.at(tmp_edge.tgt)->GetID() + "\",\n";
			ofsm << "\t\t\"weight\": \"" + utils_to_string(tmp_edge.weight) + "\"\n";
			ofsm << "\t  },\n";

			// Increment counter
			cpt++;
		}

		// Remove last comma
		pos = ofsm.tellp();
		ofsm.seekp(pos-3);
		ofsm << '\n';

		// Write edge footer
		ofsm << "\t]\n";
	}
	else
	{
		// Write node footer
		ofsm << "\t]\n";
	}

	// Write footer
	ofsm << "}";

	// Close file
	ofsm.close();

}
