/**
 *
 * BIRCH++
 * Implements BIRCH Clustering
 * 
 * DataReader.cpp
 * Data reading related functions
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-04-13
 * @version 1.0
 *
 */

#include "BIRCH.h"

/* PROTOTYPES */
void SplitContentAndStore(char *source, std::size_t nbData, vector<std::size_t> startIndices, std::size_t dimension, Label label, const char separator, double* dest, string* labels);
std::vector<Object> ArrayToObjectVector(double* dataArray, string* labelArray, std::size_t nbData, vector<std::size_t> startIndices, std::size_t dimension);


/*!
 *  \brief Load from data file
 *
 *  \param fileName  : absolute filename of the data file
 *	\param label     : position of the label if present
 *	\param separator : separator in the adta file (e.g. '\t', ';', etc.)
 *	\param nbData	 : (output) number of data-points consider
 *	\param dim	     : (output) dimension of the data-point 
 *	\param patch	 : boolean that specify if a separator is present at the end of the lines
 *
 */
std::vector<Object> LoadData(const char* filename, Label label, const char separator, std::size_t* nbData, std::size_t* dim, bool patch)
{
	// Time
	clock_t time0 = clock();

	// Load file
	ifstream file (filename, ios::in|ios::binary);
	if(!file)
	{
		std::cout << "Error: File " << filename << " not found." << std::endl;
		std::cout << "\tPress Enter to continue" << std::endl;
		std::cin.ignore();
		return std::vector<Object>();
	}

    // Get length of file
    file.seekg (0, std::ios::end);
    std::size_t length = file.tellg();
    file.seekg (0, std::ios::beg);

    // Allocate memory
    char* buffer = new char [length];

    // Read data as a block
    file.read(buffer, length);

	// Close file
	file.close();

	// Count line and get starting index
	std::size_t lines = 0;
	vector<std::size_t> lineStartIndex;
	lineStartIndex.push_back(0);

	// Parse buffer
	for(std::size_t i = 0; i < length; i++)
	{
		if (buffer[i] == '\n')
		{
			lines++;
			lineStartIndex.push_back(i+1);
		}
	}
	*nbData = lines;

	// Count dimension
	std::size_t dimension = 0;
	for(std::size_t i = 0; i < length; i++)
	{
		if (buffer[i] == separator)
		{
			dimension++;
		}
		else if (buffer[i] == '\n')
		{
			dimension++;
			break;
		}
	}

	// Patch (if last character is a separator )
	if (patch)
		dimension -= 1;
	
	// Considerer label if present
	string *labels = NULL;
	if (label != NONE)
	{
		dimension -= 1;
		labels = new string[lines];
	}
	*dim = dimension;

	// Allocate memory for data
	double *data = new double[lines * dimension];
	//double *data = (double *) malloc(lines*dimension*sizeof(double));

	// Time
	//time0 = clock() - time0;
	//std::cout << "\tFile parsed in " << (double)time0/CLOCKS_PER_SEC << " seconds" << std::endl;
	//time0 = clock();

	// Split content for futher processes
	SplitContentAndStore(buffer, lines, lineStartIndex, dimension, label, separator, data, labels);

	// Time
	//time0 = clock() - time0;
	//std::cout << "\tData stored in " << (double)time0/CLOCKS_PER_SEC << " seconds" << std::endl;
	//time0 = clock();

	// Store as Object vector
	std::vector<Object> objectVector = ArrayToObjectVector(data, labels, lines, lineStartIndex, dimension);

	// Time
	//time0 = clock() - time0;
	//std::cout << "\tObject structure created in " << (double)time0/CLOCKS_PER_SEC << " seconds" << std::endl;

	// Free memory
	delete[] buffer;
	delete[] data;
	delete[] labels;
	lineStartIndex.clear();

	return objectVector;
}

/*!
 *  \brief Convert data from string to array
 *
 *  \param source       : source string
 *	\param nbData       : number of points in the data
 *	\param startIndices : points begin position in the string
  *	\param label		: position of the label if present
 *	\param separator    : separator in the data file (e.g. '\t', ';', etc.)
 *	\param dest         : (output) array containing data
 *	\param labels		: (output) array containing label as string
 *
 */
void SplitContentAndStore(char *source, std::size_t nbData, vector<std::size_t> startIndices, std::size_t dimension, Label label, const char separator, double* dest, string* labels)
{
	char *currentSrcPtr;
	double *currentDestPtr;
	std::size_t offset = 0;
	char* currentValue;

	#pragma omp parallel for schedule(dynamic, 128) private(currentSrcPtr, currentDestPtr, offset, currentValue)
	for(__int64 lineIndex = 0; lineIndex < nbData; lineIndex++)
	{
		// Set pointer in source
		std::size_t ptr = startIndices.at(lineIndex);
		currentSrcPtr = &source[ptr];
		// Set pointer in destination
		currentDestPtr = dest + lineIndex * dimension;
		// Set offset at the begining of the line
		offset = 0;
		// Store first column if label at the first position
		if (label == -1)
		{
			currentValue = currentSrcPtr;
			while( (currentSrcPtr[offset] != separator) && (currentSrcPtr[offset] != '\n') && (currentSrcPtr[offset] != '\0') )
				offset++;
			currentSrcPtr[offset] = '\0';
			labels[lineIndex] = currentValue;
			offset++;
		}
		// Split line
		for(std::size_t k = 0; k < dimension; k++)
		{
			currentValue = currentSrcPtr + offset;
			while( (currentSrcPtr[offset] != separator) && (currentSrcPtr[offset] != '\n') && (currentSrcPtr[offset] != '\0') )
				offset++;
			currentSrcPtr[offset] = '\0';	
			currentDestPtr[k] = atof(currentValue);
			offset++;
		}	
		// Store last column if label at the last position
		if (label == 1)
		{
			currentValue = currentSrcPtr + offset;
			while( (currentSrcPtr[offset] != separator) && (currentSrcPtr[offset] != '\n') && (currentSrcPtr[offset] != '\0') )
				offset++;
			currentSrcPtr[offset] = '\0';
			labels[lineIndex] = currentValue;
			offset++;
		}
	}
	return;
}


/*!
 *  \brief Convert data from array to object vector
 *
 *  \param dataArray    : data array
 *	\param labelArray   : label array
 *	\param nbData		: number of data-points considered
 *	\param startIndices : points begin position in the string
 *	\param dimension    : separator in the data file (e.g. '\t', ';', etc.)
 *
 *	\return Data as vector of Object
 *
 */
std::vector<Object> ArrayToObjectVector(double* dataArray, string* labelArray, std::size_t nbData, vector<std::size_t> startIndices, std::size_t dimension)
{
	std::vector<Object> objectVector = std::vector<Object>(nbData);

	#pragma omp parallel for schedule(dynamic, 128) shared(dataArray, labelArray, objectVector)
	for(__int64 i=0; i<nbData; i++)
	{
		// Create an object from data
		Object newObject;
		newObject.id = i;
		newObject.tree_id = "";
		newObject.ptr = new double[dimension];
		std::copy(dataArray + i*dimension, dataArray + (i+1)*dimension, newObject.ptr);
		newObject.dimension = dimension;
		newObject.label = (labelArray == 0)? "" : labelArray[i];
		newObject.imagepath = newObject.label; 
		newObject.imagepath.erase(newObject.imagepath.find_last_not_of(" \n\r\t")+1); // right whitespace trim
		newObject.imagepath += ".jpg";
		newObject.cluster_id = -1;

		#pragma omp critical(add_object)
		{
			// Add to the vector
			objectVector[i] = newObject;
		}
	}

	return objectVector;

}

