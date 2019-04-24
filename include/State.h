#ifndef State_h
#define State_h

///
///	\class State
///
///	\brief Holder class for state information 
///
/// 	Holder class for information on nuclear states. In Coulomb 
///	excitation, each substate is treated independently, so the 
///	state only holds information common to all substates.
///

class State {

	public :
		State()	{;}
		State(double E, double J);		/*!< Construct a State of energy, E (MeV) and spin, J */
		~State() {;}
		State(const State& s);			/*!< Copy constructor */
		State& operator = (const State& s);	/*!< Assignment operator */

		void			SetPsi(double p)		{ fPsi = p;		}	/*!< Define psi, state-by-state parameter feeding into the coupling parameters*/
		void			SetEta(double e)		{ fEta = e;		}	/*!< Define eta, the wavenumber of the state  */
		void			SetPhase(int p)			{ fPhaseFactor = p;	}	/*!< Define the spin and parity phase factor */

		double			GetJ()			const	{ return fJ;		}	/*!< Return state J */
		double			GetStateEnergy()	const	{ return fStateEnergy;	}	/*!< Return state energy (MeV) */
		double			GetEta()		const	{ return fEta;		}	/*!< Return state eta */
		double			GetPsi()		const	{ return fPsi;		}	/*!< Return state psi */	
		int			GetPhase()		const	{ return fPhaseFactor;	}	/*!< Return the spin and parity phase factor */

	private :

		double			fJ;
		double			fStateEnergy;
		double			fEta;
		double			fPsi;
		int			fPhaseFactor;

};
#endif
