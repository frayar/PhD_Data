#include "main.h"
#include "export.h"
#include "irng.h"


/* PROTOTYPES */
void ProcessFile(const char* fileName, bool label, const char separator);
void SplitContentAndStore(char *source, int nbdata, vector<int> startIndices, int dimension, const char separator, double* dest);
double CalculeCPUV1(double *pfData, int iRow, int iCol, int threads, const char* outputFileName);

/* GLOBALS */
std::string g_fileName = "C:\\Fred2\\iRNG_Exact\\MIRFLICKR1M_EHD.txt";
bool g_islabel = true;
char g_separator = '\t';
ALGO g_ALGO = IRNG;
int g_NB_THREAD = 4;
EXPORT g_EXPORT_METHOD = ADJACENCIES;
std::string g_graphOutputFilename = "C:\\Fred2\\iRNG_Exact\\miflickr1m\\_irng_omp_miflickr1m.json";
std::string g_edgesOutputFilename = "C:\\Fred2\\iRNG_Exact\\miflickr1m\\_irng_omp_miflickr1m.edges";
std::string g_distancesOutputFilename = "C:\\Fred2\\iRNG_Exact\\miflickr1m\\_irng_omp_miflickr1m.distances";
std::string g_timelogFilename = "C:\\Fred2\\iRNG_Exact\\miflickr1m\\_irng_omp_miflickr1m.timelog";


/**
 *
 * MAIN 
 *
 */
int main( int argc, const char* argv[] )
{
	ProcessFile(g_fileName.c_str(), g_islabel, g_separator);
}

/**
 *
 * PROCESS ONE FILE 
 *
 */
void ProcessFile(const char* fileName, bool label, const char separator)
{
	// Load file
	ifstream file (fileName, ios::in|ios::binary);
	if(!file)
	{
		std::cout << "Error: File " << fileName << " not found." << std::endl;
		std::cout << "Press Enter to continue" << std::endl;
		std::cin.ignore();
		return;
	}

    // Get length of file
    file.seekg (0, std::ios::end);
    int length = file.tellg();
    file.seekg (0, std::ios::beg);

    // Allocate memory
    char * buffer = new char [length];

    // Read data as a block
    file.read(buffer, length);

	// Close file
	file.close();

	// Count line and get starting index
	int lines = 0;
	vector<int> lineStartIndex;
	lineStartIndex.push_back(0);

	// Parse buffer //#pragma omp parallel for schedule(dynamic, 256) shared(lines)
	for(int i = 0; i < length; i++)
	{
		if (buffer[i] == '\n')
		{
			lines++;
			lineStartIndex.push_back(i+1);
		}
	}

	// Count dimension
	int dimension = 0;
	for(int i = 0; i < length; i++)
	{
		if (buffer[i] == g_separator)
			dimension++;
		else if (buffer[i] == '\n')
		{
			dimension++;
			break;
		}
	}
	
	// Considerer last column label
	if (label)
		dimension -= 1;

	// Allocate memory for data
	double *data = new double[dimension * lines];

	// Split content for futher processes
	SplitContentAndStore(buffer, lines, lineStartIndex, dimension, separator, data);

	// Calcule RNG
	double time0 = 0.0;
	switch(g_ALGO)
	{
		case RNG:
			time0 = CalculeCPUV1(data, lines, dimension, g_NB_THREAD, g_graphOutputFilename.c_str());
			std::cout << "RNG of data(" << lines << " x " << dimension << ") with " << g_NB_THREAD << " threads processed in " << time0 << " milliseconds" << std::endl;
			break;
		case IRNG:
			time0 = Compute_iRNG(data, lines, dimension, g_NB_THREAD, g_EXPORT_METHOD, g_graphOutputFilename.c_str(), g_edgesOutputFilename.c_str(), g_distancesOutputFilename.c_str(), g_timelogFilename.c_str());
			std::cout << "iRNG of data(" << lines << " x " << dimension << ") with " << g_NB_THREAD << " threads processed in " << time0 << " seconds" << std::endl;
			break;
		default:
			break;
	}
	
	std::cout << "Press Enter to continue" << std::endl;
	std::cin.ignore();

	// Free memory
	delete[] buffer;
	delete[] data;
	lineStartIndex.clear();

	return;
}

/**
 * 
 * CONVERT FROM STRING TO ARRAY 
 *
 */
void SplitContentAndStore(char *source, int nbdata, vector<int> startIndices, int dimension, const char separator, double* dest)
{
	char *currentSrcPtr;
	double *currentDestPtr;
	int offset = 0;
	char* currentValue;

	#pragma omp parallel for schedule(dynamic, 128) private(currentSrcPtr, currentDestPtr, offset, currentValue)
	for(int lineIndex = 0; lineIndex < nbdata; lineIndex++)
	{
		// Set pointer in source
		int ptr = startIndices.at(lineIndex);
		currentSrcPtr = &source[ptr];
		// Set pointer in destination
		currentDestPtr = dest + lineIndex * dimension;
		// Set offset at the begining of the line
		offset = 0;
		// Split line
		for(int k = 0; k < dimension; k++)
		{
			currentValue = currentSrcPtr + offset;
			while( (currentSrcPtr[offset] != separator) && (currentSrcPtr[offset] != '\n') )
				offset++;
			currentSrcPtr[offset] = '\0';	
			currentDestPtr[k] = atof(currentValue);
			offset++;
		}		
	}
	return;
}


/**
 *
 * COMPUTE O(n^3) RNG 
 *
 */
double CalculeCPUV1(double *pfData, int iRow, int iCol, int threads, const char* outputFileName){

	// Attributes
	list<DoubleResultPackage> *resultsCPU = NULL;
	double time0;
	double **distance;
	double distIJ,distIK,distJK;
	double iDistance = 0.0f;
	double temp = 0.0f;
	bool isCanUse;

	// Specifies the number of threads used by default in subsequent parallel sections
	omp_set_num_threads(threads);

	// Set result structure
	DoubleResultPackage resultPackage;
	resultsCPU = new list<DoubleResultPackage>[iRow];

	// Set similarity matrix structure
	distance = new double*[iRow];
	for(int idxI = 0; idxI < iRow; ++idxI)
		distance[idxI] = new double[iRow];

	// Set clock
	time0 = clock();

	// STEP 1: Compute similarity matrix
	#pragma omp parallel for schedule(dynamic, 256) private(iDistance, temp)
	for(int idxI = 0; idxI < iRow; idxI++){
		for(int idxJ = idxI + 1; idxJ < iRow; idxJ++){	
			iDistance = 0.0f;
			for(int counter = 0, idxIi = idxI * iCol, idxJj = idxJ * iCol; counter < iCol; counter++, idxIi++, idxJj++){
				temp = pfData[idxIi] - pfData[idxJj];
				temp = temp * temp;
				iDistance += temp;
			}
			distance[idxI][idxJ] = distance[idxJ][idxI] = sqrt(iDistance);
		}

		// Set diagonal at 0
		distance[idxI][idxI] = 0.0f;
	}	
	// STEP 2: for each(i,j) isCanUse if edge(i,j) exist in RNG
	#pragma omp parallel for schedule(dynamic, 256) private(resultPackage, isCanUse,distIJ,distIK,distJK)
	for(int idxI = 0; idxI < iRow; idxI++){		
		for(int idxJ = idxI + 1; idxJ < iRow; idxJ++){
			isCanUse = true;
			distIJ = distance[idxI][idxJ];	
			for(int idxK = 0; idxK < iRow; idxK++){
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
				resultPackage.Distance = distIJ;
				resultPackage.Idx  = idxJ;
				resultsCPU[idxI].push_back(resultPackage);	
			}
		}
	}	

	// Get RNG time
	time0 = clock() - time0;
	//time0 = time0/(double)1000.0;

	// Export graph
	ExportGraph(resultsCPU, outputFileName, iRow, g_EXPORT_METHOD);

	// Debug: nb edges
	
	size_t nbEdges = 0;
	for (int i=0; i<iRow; i++)
		nbEdges += resultsCPU[i].size();
	std::cout << "Total edges = " << nbEdges << std::endl;
	

	// Free memory
	for(int idxI = 0; idxI < iRow; ++idxI)
		delete[] distance[idxI];
	delete[] distance;
	delete[] resultsCPU;
	return time0;
}
