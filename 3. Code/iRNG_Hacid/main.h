#ifndef MAIN_H
#define MAIN_H

	#include <iostream>
	#include <time.h>
	#include <omp.h>
	#include <list>
	#include <map>
	#include <vector>
	#include <fstream>
	#include <math.h> 
	#include <omp.h>
	#include <limits>		// std::numeric_limits<float>::infinity()
	#include <algorithm>    // std::sort

	using namespace std;

	typedef struct DoubleResultPackage{
		int Idx;
		double Distance;

		// WARNING: NOT MULTIGRAPH !!!
		bool operator == (const DoubleResultPackage &other)
		{
			return (other.Idx == Idx);
		}

		void print()
		{
			std::cout << "ID = " << Idx << " - distance = " <<  Distance << std::endl;
		}

	} DoubleResultPackage;

	enum ALGO {RNG, IRNG};
	enum EXPORT { ADJACENCIES, JSON, GEXF, TLP };
	enum UPDATE { O_N_CUBE, TWO_PASS };

	const int IRNG_EPSILON = 1;

#endif

