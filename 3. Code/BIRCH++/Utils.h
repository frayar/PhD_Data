/**
 *
 * BIRCH++ 
 * Implements BIRCH Clustering
 *
 * utils functions
 *
 *
 * @author	Frédéric RAYAR <frederic.rayar@univ-tours.fr>
 * @date	2015-05-29
 * @version 1.0
 *
 */

#ifndef __UTILS_H__
#define	__UTILS_H__


	#include <sstream>

	namespace 
	{

		/*!
		 *
		 *	\brief Convert from std::size_t to std::string
		 *
		 *	\param value : std::size_t to convert
		 *
		 *	\returns std::string of the std::size_t
		 *
		 */
		std::string utils_to_string(std::size_t value)
		{
			std::stringstream ss;
			ss << value;

			return ss.str();
		}

		/*!
		 *
		 *	\brief Convert from int to std::string
		 *
		 *	\param value : int to convert
		 *
		 *	\returns std::string of the int
		 *
		 */
		std::string utils_to_string(int value)
		{
			std::stringstream ss;
			ss << value;

			return ss.str();
		}


		/*!
		 *
		 *	\brief Convert from double to std::string
		 *
		 *	\param value : double to convert
		 *
		 *	\returns std::string of the std::size
		 *
		 */
		std::string utils_to_string(double value)
		{
			std::stringstream ss;
			ss << value;

			return ss.str();
		}

	}



#endif