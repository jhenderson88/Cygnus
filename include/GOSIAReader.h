#ifndef GOSIAReader_h
#define GOSIAReader_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Nucleus.h"
#include "ExperimentalInput.h"


class GOSIAReader {

	public:

		GOSIAReader(Nucleus*,const char*);	
		~GOSIAReader()				{;}

		void	Clear()				{ gosiaData.clear();		}
		void	ReadGOSIAFile(const char*);			/*!< Read formatted GOSIA file and extract yields */

		std::vector<ExperimentData>		GetGOSIAData()		{ return gosiaData;	}	/*!< Return vector of ExperimentData containing GOSIA yields*/

		void	Print()	const;

	private:

		Nucleus					fNucleus;	/*!< Nucleus object used to define level assignments */
		std::vector<ExperimentData>		gosiaData;	/*!< Vector of ExperimentData objects to hold the data that is read in from the file */


};

#endif
