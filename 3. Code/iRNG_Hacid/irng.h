#ifndef IRNG_H
#define IRNG_H

	/* PROTOTYPES */
	double Compute_iRNG(double *pfData, int iRow, int iCol, int threads, EXPORT exportMethod, const char* graphOutputFilename, const char* edgesOutputFilename, const char* distancesOutputFilename, const char* timelogFilename);
	double insert_point(__int64 pointID, double *data, __int64 nbdata, __int64 dimension, int threads, __int64 *nearest, __int64 *farthest, __int64 *sr_candidates, std::map<__int64, double> *resultsCPU);
	__int64 GetNearest(std::map<__int64, double> neighbors);
	__int64 GetFarthest(std::map<__int64, double> neighbors);
	double my_distance(double *data, __int64 dimension, __int64 p1_id, __int64 p2_id);

	// load_raw_points
	// load_similarity_matrix
	// insert_points
	// insert_point
	// remove_point
	// insert_point_optimized (dmin/dmax)
	// remove_point_optimized (heuristic)

#endif