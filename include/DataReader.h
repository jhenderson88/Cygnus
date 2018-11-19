#ifndef DataReader_h
#define DataReader_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Nucleus.h"
#include "ExperimentalInput.h"

///
///	\class DataReader
///
///	\brief Simple class used to read formatted experimental yields
///	and put it in a vector of ExperimentData objects for fitting.	
///

class DataReader {

	public:
		DataReader(Nucleus*,const char*);			/*!< Constructor with Nucleus (for level assignments) and filename of formatted data */
		~DataReader()				{;	}

		void	ReadDataFile(const char*);			/*!< Read formatted yield file and add data for ExperimentData vector */

		std::vector<ExperimentData>		GetExperimentData()		{ return exptData;	}	/*!< Return vector of ExperimentData */

	private:
		Nucleus					fNucleus;	/*!< Nucleus object used to define level assignments */
		std::vector<ExperimentData>		exptData;	/*!< Vector of ExperimentData objects to hold the data that is read in from the file */

};

#endif
