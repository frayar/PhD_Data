#include "main.h"
#include "export.h"


/*
 * EXPORT RNG - FOR RNG ALGO
 * TODO: Change resultsCPU structure of RNG algorithm to std::map
 *        and remove this duplicate
 *
 */
void ExportGraph(list<DoubleResultPackage>* pList, const char *filePath, int row, EXPORT method)
{
	if (method == GEXF)
	{
	
		ofstream ofsm(filePath , ios_base::trunc);

		// Header
		ofsm << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
		ofsm << "<gexf xmlns=\"http://www.gexf.net/1.2draft\" version=\"1.2\">\n";
		ofsm << "\t<graph mode=\"static\" defaultedgetype=\"directed\">\n";
		ofsm << "\t\t<nodes>\n";

		for(__int64 idxI = 0; idxI < row; ++idxI)
			ofsm <<  "\t\t\t<node id=\"" << idxI << "\" label=\"" << idxI << "\" />\n";
		
		ofsm << "\t\t</nodes>\n";
		ofsm << "\t\t<edges>\n";

		// Edges
		int cpt = 0;
		for(__int64 idxI = 0; idxI < row; ++idxI)
		{
			for(list<DoubleResultPackage>::iterator it = pList[idxI].begin(); it != pList[idxI].end(); it++)
			{
				if ((*it).Idx > idxI)
				{
					ofsm << "\t\t\t<edge id=\"" << cpt << "\" source=\"" << idxI << "\" target=\"" << (*it).Idx << "\" weight=\"" << (*it).Distance << "\"  />\n";
					cpt++;
				}
			}
		}

		ofsm << "\t\t</edges>\n";
		ofsm << "\t</graph>\n";
		ofsm << "</gexf>\n";
		
		ofsm.close();

	}
	else
	{
		//Build edges list
		std::multimap<__int64, __int64> edges;
		for(int idxI = 0; idxI < row; ++idxI)
		{
			for(list<DoubleResultPackage>::iterator it = pList[idxI].begin(); it != pList[idxI].end(); it++)
				if ((*it).Idx > idxI)
					edges.insert(std::pair<__int64, __int64>(idxI,(*it).Idx));
				else
					edges.insert(std::pair<__int64, __int64>((*it).Idx,idxI));
		}

		// Export edges list
		ofstream ofsm(filePath , ios_base::trunc);
		for(std::multimap<__int64, __int64>::iterator it = edges.begin(); it != edges.end(); it++)
			ofsm << it->first << '\t' << it->second << endl;
	}


}




/**
 *
 * EXPORT GRAPH
 *
 */
void ExportGraph(std::map<__int64, double>* pList, const char *filePath, int row, EXPORT method){
	switch (method)
	{
		case ADJACENCIES:
			WriteAdjacencies(pList, filePath, row);
			break;
		case JSON:
			WriteJSONResults(pList, filePath, row);
			break;
		case GEXF:
			WriteGEXFResults(pList, filePath, row);
			break;
		case TLP:
			WriteTLPResults(pList, filePath, row);
			break;
		default:
			WriteAdjacencies(pList, filePath, row);
			break;
	}
}




/**
 *
 * EXPORT RNG ADJACENCIES (for NetworkX)
 *
 */
void WriteAdjacencies(std::map<__int64, double>* pList, const char *filePath, int row){
	ofstream ofsm(filePath , ios_base::trunc);
	for(__int64 idxI = 0; idxI < row; ++idxI)
	{
		for(std::map<__int64, double>::iterator it = pList[idxI].begin(); it != pList[idxI].end(); it++)
			if (it->first > idxI)
				ofsm << idxI << '\t' << it->first << endl;
	}
	ofsm.close();
}




/**
 *
 * EXPORT RNG JSON for oriented graph
 * Handle case when for a given edge between i and j, the structure pList[i] contains j AND pList[j] contains i 
 *
 */
void WriteJSONResults(std::map<__int64, double>* pList, const char *filePath, int row)
{

	ofstream ofsm(filePath , ios_base::trunc);
	ofsm << "{\n   \"directed\": false,\n   \"graph\": [],\n   \"nodes\": [\n";

	// Nodes
	for(__int64 idxI = 0; idxI < row-1; ++idxI)
	{
		ofsm << "      {\n         \"id\": " << idxI << "\n      },\n";
	}
	// Last node (no comma)
	ofsm << "      {\n         \"id\": " << row-1 << "\n      }\n";
	ofsm << "   ],\n   \"links\": [\n";

	// Edges
	for(int idxI = 0; idxI < row-1; ++idxI)
	{
		for(std::map<__int64, double>::iterator it = pList[idxI].begin(); it != pList[idxI].end(); it++)
			if (it->first > idxI)
				ofsm << "      {\n         \"source\": " << idxI << ",\n         \"target\": " << it->first << ",\n         \"weight\": " << it->second << "\n      },\n";
	}

	// Last edges (no comma)
	// TODO: Find/Replace last occurence of "}," by "}"

	ofsm << "   ],\n   \"multigraph\": false\n}";
	ofsm.close();
}




/**
 *
 * EXPORT RNG GEXF for TULIP/GEPHI
 *
 */
void WriteGEXFResults(std::map<__int64, double>* pList, const char *filePath, int row)
{
	ofstream ofsm(filePath , ios_base::trunc);

	// Header
	ofsm << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	ofsm << "<gexf xmlns=\"http://www.gexf.net/1.2draft\" version=\"1.2\">\n";
	ofsm << "\t<graph mode=\"static\" defaultedgetype=\"directed\">\n";
	ofsm << "\t\t<nodes>\n";

	for(__int64 idxI = 0; idxI < row; ++idxI)
		ofsm <<  "\t\t\t<node id=\"" << idxI << "\" label=\"" << idxI << "\" />\n";
	
	ofsm << "\t\t</nodes>\n";
	ofsm << "\t\t<edges>\n";

	// Edges
	int cpt = 0;
	for(__int64 idxI = 0; idxI < row; ++idxI)
	{
		for(std::map<__int64, double>::iterator it = pList[idxI].begin(); it != pList[idxI].end(); it++)
		{
			if (it->first > idxI)
			{
				ofsm << "\t\t\t<edge id=\"" << cpt << "\" source=\"" << idxI << "\" target=\"" << it->first << "\" weight=\"" << it->second << "\"  />\n";
				cpt++;
			}
		}
	}

	ofsm << "\t\t</edges>\n";
	ofsm << "\t</graph>\n";
	ofsm << "</gexf>\n";
	
	ofsm.close();
}




/**
 *
 * EXPORT RNG TLP for TULIP
 *
 */
void WriteTLPResults(std::map<__int64, double>* pList, const char *filePath, int row)
{
	ofstream ofsm(filePath , ios_base::trunc);

	// Header
	ofsm << "(tlp \"2.0\"\n";
	ofsm << "(date \"10-23-2006\")\n";
	ofsm << "(author \"frederic.rayar\")\n";

	// Nodes
	ofsm << ";(nodes <node_id> <node_id> ...)\n";
	ofsm << "(nodes";
	for(__int64 idxI = 0; idxI < row; ++idxI)
		ofsm <<  " " << idxI;
	ofsm << ")\n";

	// Edges
	ofsm << ";(edge <edge_id> <source_id> <target_id>)\n";
	__int64 cpt = 0;
	for(__int64 idxI = 0; idxI < row; ++idxI)
	{
		for(std::map<__int64, double>::iterator it = pList[idxI].begin(); it != pList[idxI].end(); it++)
		{
			if (it->first > idxI)
			{
				ofsm << "(edge " << cpt << " " << idxI << " " << it->first << ")\n";
				cpt++;
			}
		}
	}

	// Weights
	ofsm << "\n";
	ofsm << "(property  0 double \"viewMetric\"\n";
	ofsm << "(default \"0\" \"0\")\n";
	cpt = 0;
	for(__int64 idxI = 0; idxI < row; ++idxI)
	{
		for(std::map<__int64, double>::iterator it = pList[idxI].begin(); it != pList[idxI].end(); it++)
		{
			if (it->first > idxI)
			{
				ofsm << "(edge " << cpt << " \"" << it->second << "\")\n";
				cpt++;
			}
		}
	}

	// Footer
	ofsm << ")";
	ofsm.close();
}