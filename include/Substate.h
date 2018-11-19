#ifndef Substate_h
#define Substate_h

#include <cstddef>
#include <vector>
#include "Connection.h"

///
///	\class Substate
///
///	\brief Holder class for substate information
///
/// 	Holder class for information on magnetic substates including M, 
///	connection information, etc.                 
///
///	Excitation calculations are performed substate-by-substate, so
///	this class (rather than the State class) contains the majority
///	of the important coupling information.
///
///	Contained within this class is a list of connections to other 
///	substates, defined by the Connections class.
///
///	Also included is an index for the State object, of which the
///	substate is a... substate.
///

class Substate {

	public :
		Substate() {;}
		/// 
		///	Constructor defining the magnetic quantum number, M, the substate index, ss, and the 
		///	parent state index, s.
		///
		Substate(double M, int ss, int s);
		~Substate() {;}
		Substate(const Substate& s);			/*!< Copy constructor */
		Substate& operator = (const Substate& s);	/*!< Assignment operator */

		int			GetStateIndex()		const	{ return fStateIndex;			}	/*!< Return the parent state index */
		void			SetStateIndex(int i)		{ fStateIndex = i;			}	/*!< Define the parent state index */

		double			GetM()			const	{ return fM;				}	/*!< Return the magnetic quantum number */
		void			SetM(double m)			{ fM = m;				}	/*!< Define the magnetic quantum number */

		unsigned int		GetNconnections()	const	{ return fConnections.size();		}	/*!< Return the number of connections to other substates */
		Connection		GetConnection(int i)	const	{ return fConnections.at(i);		}	/*!< Return the connection (indexed i) to another substate */
	
		void			AddConnection(Connection c)	{ fConnections.push_back(c);		}	/*!< Define a new connection to another substate */

	private:
		int 			fSubstateIndex;
		int 			fStateIndex;
		double			fM;
		std::vector<Connection>	fConnections;

};

#endif
