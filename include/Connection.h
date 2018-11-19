#ifndef Connection_h
#define Connection_h

#include <cstddef>
#include <vector>

///
///	\class Connection
///
///	\brief 	Holder for connections between substates (see Substate class). 
///	Contains information relevant to point Coulomb excitation calculations.
///	
///	Contains vectors of connection information. Only one Connection is defined 
///	for any connection between substates, which contains all multipolarities in 
///	that connection. 
///
///	Size of the vectors is defined using SetMaxLambda(), for convenience the 
///	vectors between all substates are given the same length.
///
///	Because each connection is attached to a Substate, the fConnectedState 
///	index defines the other end of the connection only (i.e. the Connection
///	does not contain information about the index of the state it is attached to).
///

class Connection {

	public:
		Connection() {;}	
		~Connection() {;}
		Connection(const Connection& c);		/*!< Copy constructor */
		Connection& operator = (const Connection& c);	/*!< Assignment operator */

		void	SetConnectedState(int s)	{ fConnectedState = s;			}	/*!< Define the index of the substate the holder substate is connected to by this object */

		int	GetConnectedState()		{ return fConnectedState;		}	/*!< Return the index of the substate the holder substate is connected to by this object */
		bool	IsSet(int i)			{ return fLambda.at(i);			}	/*!< Returns true if the connection for the multipolarity defined by i is being used, else false */
		double	GetPsi(int i)			{ return fPsi.at(i);			}	/*!< Returns the psi value of the connection */
		double	GetXi(int i)			{ return fXi.at(i);			}	/*!< Returns the xi value of the connection */
		double	GetZeta(int i)			{ return fZeta.at(i);			}	/*!< Returns the zeta value of the connection */
		int	GetN()			const 	{ return fXi.size();			}	/*!< Returns the number of multipolarites defined (including empty) */

	
		int	GetMaxLambda()		const	{ if(fPsi.size() == fXi.size() && fPsi.size() == fZeta.size()) 
								return fPsi.size();			
							  else
								return 0;
							}						/*!< Safer check on the number of multipolarities defined */

		void	SetMaxLambda(int L)		{ 
								fPsi.clear();
								fXi.clear();
								fZeta.clear();
								fLambda.clear();
								fPsi.resize((unsigned int)L); 
								fXi.resize((unsigned int)L); 
								fZeta.resize((unsigned int)L);	
								fLambda.resize((unsigned int)L);
								for(int i=0;i<L;i++)
									fLambda.at(i) = false;
							}						/*!< Define maximum multipolarity: sets the size of the vectors in the class */

		void	AddLambda(int l, double xi, double psi, double zeta);				/*!< Define xi, psi and zeta for multipolarity defined by l */
		void	Clear();									/*!< Clears all vectors */

	private:
		std::vector<bool>	fLambda;							/*!< Vector of bools to check whether the mulitpolarity is defined */
		std::vector<double>	fPsi;								/*!< Vector of Psi values */
		std::vector<double>	fXi;								/*!< Vector of Xi values */
		std::vector<double>	fZeta;								/*!< Vector of Zeta values */
		int			fConnectedState;						/*!< Index of the substate the holder substate is connected to through this object */

};
#endif 
