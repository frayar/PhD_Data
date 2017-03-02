/**
 *
 * BIRCH++
 * 
 * Implements BIRCH Clustering
 * Original paper: "BIRCH: an efficient data clustering method for very large databases", Zhang T., Ramakrishnan R. and Livny M. (1996)
 * Inspired by birch-clustering-algorithm by Taesik Yoon (otterrrr@gmail.com) (https://code.google.com/p/birch-clustering-algorithm/)
 *
 *
 * Usage : birch <data-points-file>
 *           arg1 - data-points file
 *
 *
 * The program consists in 5 mains steps:
 * - Load data-points
 * - Build a CFTree from data-points
 * - Compact the CFTree
 * - Cluster the leaves entries
 * - Assign cluster to data-points
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-04-13
 * @version 1.0
 *
 */

#include "BIRCH.h"
#include "CFTree.h"


/* PROTOTYPES */
std::vector<Object> LoadData(const char* filename, Label label, const char separator, std::size_t* nbData, std::size_t* dim, bool patch);
void ExportTree_JSON(CFTree *tree, std::string filename);
void ExportLeaves_JSON(CFTree *tree, std::string filename);
void ExportGraph_JSON(std::string level, std::vector<CFEntry*> nodes, std::vector<Edge> edges, std::string filename);
void ComputeRepresentatives(CFTree* tree, double time0);
void GetHierarchicalRNG(CFTree* tree, double time0);
void GetMultilevelRNG(CFTree* tree, double time0);
void ExportHierarchicalRNG(CFTree* tree, double time0);
void ExportLeaves(CFTree* tree, double time0);
void ExportMultilevelRNG(CFTree* tree, double time0);


/* GLOBALS */
// Wang
/*
std::string g_filename		= "C:\\Users\\rayar\\Desktop\\2k15\\GraphViewer_v1\\data\\Wang\\Wang_CLD.txt";
Label g_label				= END;
char g_separator			= '\t';
bool g_patch				= false;
std::string g_outputfile	= "C:\\Users\\rayar\\Desktop\\2k15\\GraphViewer_v1\\data\\Wang\\Wang_CLD_HRNG.json";
*/

// NGA Images - EHD
/*
std::string g_filename		= "C:\\Users\\rayar\\Desktop\\2k15\\GraphViewer_v1\\data\\NGA_Images\\NGA_Images_EHD.txt";
Label g_label				= START;
char g_separator			= '\t';
bool g_patch				= false;
std::string g_outputfile	= "C:\\Users\\rayar\\Desktop\\2k15\\GraphViewer_v1\\data\\NGA_Images\\NGA_Images_EHD_HRNG.json";
*/

// NGA Images - CLD
/*
std::string g_filename		= "C:\\Users\\rayar\\Desktop\\2k15\\GraphViewer_v1\\data\\NGA_Images\\NGA_Images_CLD.txt";
Label g_label				= START;
char g_separator			= ';';
bool g_patch				= true;
std::string g_outputfile	= "C:\\Users\\rayar\\Desktop\\2k15\\GraphViewer_v1\\data\\NGA_Images\\NGA_Images_CLD_HRNG.json";
*/

// NASA APOD Images - CLD
/*
std::string g_filename		= "D:\\Dataset\\NASA-APOD\\NASA_APOD_resized_CLD.txt";
Label g_label				= START;
char g_separator			= ';';
bool g_patch				= true;
std::string g_outputfile	= "D:\\Dataset\\NASA-APOD\\NASA_APOD_resized_CLD_HRNG.json";
*/

// MIRFLICKR-25000 - CLD
/*
std::string g_filename		= "D:\\Dataset\\MIRFLICKR-25000\\MIRFLICKR-25000_CLD.txt";
Label g_label				= START;
char g_separator			= ';';
bool g_patch				= true;
std::string g_outputfile	= "D:\\Dataset\\MIRFLICKR-25000\\MIRFLICKR-25000_CLD_HRNG.json";
*/

// MIRFLICKR-25000_5k - CLD
/*
std::string g_filename		= "D:\\Dataset\\MIRFLICKR-25000\\MIRFLICKR-25000_CLD_5k.txt";
Label g_label				= START;
char g_separator			= ';';
bool g_patch				= true;
std::string g_outputfile	= "D:\\Dataset\\MIRFLICKR-25000\\MIRFLICKR-25000_CLD_5k_HRNG.json";
*/

// MIRFLICKR-1M
/*
std::string g_filename			= "D:\\Dataset\\MIRFLICKR1M\\Features\\Edge Histogram Descriptor\\features_edgehistogram\\MIRFLICKR1M_EHD.txt";
Label g_label				= NONE;
char g_separator			= '\t';
bool g_patch				= true;
std::string g_outputfile	= "D:\\Dataset\\MIRFLICKR1M\\Features\\Edge Histogram Descriptor\\features_edgehistogram\\MIRFLICKR1M_EHD_MLRNG.json";
*/

// MIRFLICKR-1M - Unix B
/*
std::string g_filename		= "C:\\temp\\fred\\MIRFLICKR1M_EHD.txt";
Label g_label				= NONE;
char g_separator			= '\t';
bool g_patch				= true;
std::string g_outputfile	= "C:\\temp\\fred\\MIRFLICKR1M_EHD_MLRNG.json";
*/

// MNIST
/*
//std::string g_filename		= "C:\\Users\\rayar\\Desktop\\2k15\\MNIST\\MNIST_DT_pruning_OPT_INTER_CLUST_order_4_MEAN_cg_points.txt";
//Label g_label				= END;
//char g_separator			= '\t';
//bool g_patch				= false;
//std::string g_outputfile	= "C:\\Users\\rayar\\Desktop\\2k15\\MNIST\\MNIST_DT_pruning_OPT_INTER_CLUST_order_4_MEAN_cg_points__MLRNG.json";
*/

// Castle
/*
std::string g_filename		= "C:\\Users\\rayar\\Desktop\\2k15\\CESR_Castle\\CESR_Castle_EHD_dot_label.txt";
Label g_label				= END;
char g_separator			= '\t';
bool g_patch				= false;
std::string g_outputfile	= "C:\\Users\\rayar\\Desktop\\2k15\\CESR_Castle\\CESR_Castle_EHD_dot_label_MLRNG.json";
*/

// 2016 - DEBUG - WANG
/*
std::string g_filename		= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\Wang_CLD.txt";
Label g_label				= END;
char g_separator			= '\t';
bool g_patch				= false;
std::string g_outputfile	= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\Wang_CLD_HRNG.json";
std::string g_ML_outputfile	= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\Wang_CLD_MLRNG.json";
std::string g_l_outputfile = "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\Wang_CLD_Leaves.json";
*/
// 2016 - DEBUG - NASA APOD
/*
std::string g_filename		= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\NASA_APOD_resized_CLD.txt";
Label g_label				= START;
char g_separator			= ';';
bool g_patch				= true;
std::string g_outputfile	= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\NASA_APOD_resized_CLD_HRNG.json";
std::string g_ML_outputfile	= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\NASA_APOD_resized_CLD_MLRNG.json";
std::string g_l_outputfile = "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\NASA_APOD_resized_CLD_Leaves.json";
*/

// 2016 - ImageNet
/*
std::string g_filelist		= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\ImageNet\\synsets.txt";
std::string g_filedir		= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\ImageNet\\synsets\\";
std::string g_synsetsTime	= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\ImageNet\\synsets.times";
std::string g_imagesTime	= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\ImageNet\\imagenet.times";
std::string g_backupfile	= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\ImageNet\\ImageNet_CLD_HRNG_";
std::string g_outputfile	= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\ImageNet\\ImageNet_CLD_HRNG.json";
std::string g_ML_outputfile	= "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\ImageNet\\ImageNet_CLD_MLRNG.json";
std::string g_l_outputfile  = "C:\\Users\\rayar\\Desktop\\2k16\\BIRCH++\\Debug\\ImageNet\\ImageNet_CLD_Leaves.json";
Label g_label				= START;
char g_separator			= ';';
bool g_patch				= true;
std::size_t	dimension		= 192;
*/

// 2016 - ImageNet Unix
std::string g_filelist		= "/home/fred/Documents/ImageNet/synsets.txt";
std::string g_filedir		= "/home/fred/Documents/ImageNet_MPEG-7_CLD_whole/";
std::string g_synsetsTime	= "/home/fred/Documents/ImageNet/synsets.times";
std::string g_imagesTime	= "/home/fred/Documents/ImageNet/imagenet.times";
std::string g_backupfile	= "/home/fred/Documents/ImageNet/ImageNet_CLD_HRNG_";
std::string g_outputfile	= "/home/fred/Documents/ImageNet/ImageNet_CLD_HRNG.json";
std::string g_ML_outputfile	= "/home/fred/Documents/ImageNet/ImageNet_CLD_MLRNG.json";
std::string g_l_outputfile  = "/home/fred/Documents/ImageNet/ImageNet_CLD_Leaves.json";
Label g_label				= START;
char g_separator			= ';';
bool g_patch				= true;
std::size_t	dimension		= 192;


bool g_DEBUG				= true;


/**
 * Process one synset
 */
std::size_t ProcessSynset(std::string filename, CFTree *tree)
{
	// Variables
	std::size_t nbData = 0;
	std::ofstream ofsm;
	struct timeval time_start;
	struct timeval time_stop; 
	long time_seconds;
	long time_milliseconds;
	long time_microseconds;

	CFEntry entry;

	// Load data
	vector<Object> objects = LoadData(filename.c_str(), g_label, g_separator, &nbData, &dimension, g_patch);

	// Open imagenet times files
	ofsm.open (g_imagesTime, std::ofstream::out | std::ofstream::app);

	// Incrementally insert all data-point in the CFTree
	for (std::size_t i = 0; i < nbData; i++)
	{
		// Time
		gettimeofday(&time_start, NULL);

		// Create an CFEntry from the data-point
		entry = CFEntry(&objects.at(i));

		// Insert the CFEntry in the tree
		tree->InsertAndUpdate(entry);

		// Time
		gettimeofday(&time_stop, NULL);
		time_seconds  = time_stop.tv_sec - time_start.tv_sec;
		time_microseconds = time_stop.tv_usec - time_start.tv_usec;
		time_milliseconds = ((time_seconds) * 1000 + time_microseconds/1000.0) + 0.5;
		ofsm << time_milliseconds << "\n";
	}

	// Free memory
	for(std::vector<Object>::iterator it = objects.begin(); it != objects.end(); it++)
		delete[] it->ptr;
	objects.clear();

	// Close file
	ofsm.close();

	return nbData;
}

/**
 * MAIN
 * Process BIRCH Clustering
 */
int main(int argc, char* argv[])
{
	// Variables
	std::vector<std::string> synsets;
	std::string line;
	std::string synset;
	std::size_t cpt = 0;
	std::size_t nb_objects = 0;
	std::string backupfile;
	std::ofstream ofsm;
	struct timeval time_start_0;
	struct timeval time_stop_0; 
	long time_seconds_0;
	struct timeval time_start_1;
	struct timeval time_stop_1; 
	long time_seconds_1;

	//_crtBreakAlloc = 914;

	// Open synset list
	ifstream file (g_filelist, ios::in|ios::binary);
	if(!file)
	{
		std::cout << "Error: The file " << g_filelist.c_str() << " not found." << std::endl;
		std::cout << "\tPress Enter to continue" << std::endl;
		std::cin.ignore();
		return -1;
	}

	// Read and store all synset ( Note: line.length()-1 if we have to remove the '\r' windows character)
	while (std::getline(file, line, '\n'))
		synsets.push_back(line.erase(line.length()));

	// Close file
	file.close();

	// Time
	gettimeofday(&time_start_0, NULL);

	// Initialise CF tree
	CFTree tree = CFTree(dimension);

	// loop on each synset file
	for(std::vector<std::string>::iterator it = synsets.begin(); it != synsets.end(); it++)
	{
		// Increase synset counter
		cpt++;

		// User notification 
		std::cout << cpt << std::endl;

		// Get current synset filename
		synset = *it;

		// Open synsets times file
		ofsm.open (g_synsetsTime, std::ofstream::out | std::ofstream::app);

		// Time
		gettimeofday(&time_start_1, NULL);

		// Process the synset file
		nb_objects = ProcessSynset(g_filedir + synset, &tree);

		// Time
		gettimeofday(&time_stop_1, NULL);
		time_seconds_1  = time_stop_1.tv_sec - time_start_1.tv_sec;

		// Export computation time
		ofsm << synset << "\t" << nb_objects << "\t" << time_seconds_1 << "\n";

		// Close time file
		ofsm.close();

		// Scheduled backup
		/*
		if (cpt % 250 == 0)
		{
			// Build filename
			backupfile = g_backupfile + utils_to_string(cpt) + ".json";

			// Export backup HRNG
			ExportTree_JSON(&tree,backupfile);
		}
		*/
	}

	// Time
	gettimeofday(&time_stop_0, NULL);
	time_seconds_0  = time_stop_0.tv_sec - time_start_0.tv_sec;
	std::cout << "\nProcessed in " << time_seconds_0  << " seconds." << std::endl;

	// Export HRNG
	double time0 = clock();
	ExportHierarchicalRNG(&tree, time0);
	/*
	// Export Leaves
	time0 = clock();
	ExportLeaves(&tree, time0);

	// Build MLRNG
	time0 = clock();
	GetMultilevelRNG(&tree, time0);

	// Export HRNG
	time0 = clock();
	ExportMultilevelRNG(&tree, time0);
	*/

	// Notification
	std::cout << "\nDone!" << std::endl;
	std::cout << "\nPress 'Enter' to continue" << std::endl;
	std::cin.ignore();

	// Memory leak detection
	//_CrtDumpMemoryLeaks();

	return 0;
}


/**
 * Batch - Deprecated
 * Compute representatives for each CFEntries
 */
void ComputeRepresentatives(CFTree* tree, double time0)
{
	// Notification
	std::cout << "\nComputing CFEntry representatives..." << std::endl;

	// Create the Multilevel RNG
	tree->ComputeRepresentatives();

	// Time
	time0 = clock() - time0;
	std::cout << "\tRepresentatives computed in " << (double)time0/CLOCKS_PER_SEC << " seconds" << std::endl;
}


/**
 * Batch - Deprecated
 * Compute Hierachical RNG
 * Here, a hierachical RNG correspond to the RNG yield for each node entries
 */
void GetHierarchicalRNG(CFTree* tree, double time0)
{
	// Notification
	std::cout << "\nCreating Hierarchical RNG ..." << std::endl;

	// Create the Multilevel RNG
	tree->CreateHierarchicalRNG();
}


/**
 * Export Hierachical RNG
 * Here, a hierachical RNG correspond to the RNG yield for each node entries
 */
void ExportHierarchicalRNG(CFTree* tree, double time0)
{
	// Notification
	std::cout << "\nExporting Hiearchical RNG to JSON ..." << std::endl;

	// Export the ML-RNG
	ExportTree_JSON(tree,g_outputfile); 

	// Time
	time0 = clock() - time0;
	std::cout << "\tHiearchical RNG JSON exported in " << (double)time0/CLOCKS_PER_SEC << " seconds" << std::endl;
}

/**
 * Export Leaves
 * They might not be at the same level.
 */
void ExportLeaves(CFTree* tree, double time0)
{
	// Notification
	std::cout << "\nExporting Hiearchical RNG to JSON ..." << std::endl;

	// Export the ML-RNG
	ExportLeaves_JSON(tree,g_l_outputfile); 

	// Time
	time0 = clock() - time0;
	std::cout << "\tHiearchical RNG JSON exported in " << (double)time0/CLOCKS_PER_SEC << " seconds" << std::endl;
}


/**
 * Batch - Compute Multilevel RNG
 * Here, a multilevel RNG correspond to the RNG yield for each level node entries 
 */
void GetMultilevelRNG(CFTree* tree, double time0)
{
	// Notification
	std::cout << "\nCreating Multilevel RNG ..." << std::endl;

	// Compute tree height
	std::size_t height = tree->ComputeHeight();

	// Notification
	std::cout << "\tTree height = " << height << std::endl;

	// Create the Multilevel RNG
	tree->CreateMultilevelRNG();

	// Time
	time0 = clock() - time0;
	std::cout << "\tMultilevel RNG created in " << (double)time0/CLOCKS_PER_SEC << " seconds" << std::endl;
}


/**
 * Export Multilevel RNG
 * Here, a multilevel RNG correspond to the RNG yield for each level node entries 
 */
void ExportMultilevelRNG(CFTree* tree, double time0)
{
	// Notification
	std::cout << "\nExporting Multilevel RNG to JSON ..." << std::endl;

	// Export the ML-RNG
	size_t pos = g_ML_outputfile.find_last_of('.');
	std::string tmp_outputfile;
	for (std::size_t l=0; l<tree->GetHeight() + 1 ; l++)
	{
		tmp_outputfile = g_ML_outputfile;
		tmp_outputfile.insert(pos, "_" + utils_to_string(l));
		ExportGraph_JSON(utils_to_string(l), tree->GetRngsNodes()->at(l), tree->GetRngsEdges()->at(l), tmp_outputfile);
	}

	// Time
	time0 = clock() - time0;
	std::cout << "\tMultilevel RNG JSON exported in " << (double)time0/CLOCKS_PER_SEC << " seconds" << std::endl;
}


