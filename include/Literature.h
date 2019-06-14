#ifndef Literature_h
#define Literature_h

///
///	\class LitLifetime
///
///	\brief Holder for literature lifetime information for fitting
///

class LitLifetime{			// Literature state lifetime information holder

	public:
		LitLifetime()	{;}
		///
		///	Construct literature lifetime with symmetric uncertainties:\n
		///	State index = i \n
		///	State lifetime = l \n
		///	Uncertainty = e 
		///
		LitLifetime(int i, double l, double e)			{ StateIndex = i; Lifetime = l; UpUncertainty =  e; DnUncertainty =  e; 	}
		///
		///	Construct literature lifetime with asymmetric uncertainties:\n
		///	State index = i \n
		///	State lifetime = l \n
		///	Positive uncertainty = ue 
		///	Negative uncertainty = de
		///
		LitLifetime(int i, double l, double ue, double de)	{ StateIndex = i; Lifetime = l; UpUncertainty = ue; DnUncertainty = de; 	}		// State index, lifetime, down error, up error
		~LitLifetime()	{;}
		LitLifetime(const LitLifetime& lt);					/*!< Copy constructor */
		LitLifetime& operator = (const LitLifetime& lt);			/*!< Assignment operator */

		double		GetIndex()		const		{ return StateIndex;		}	/*!< Return state index */
		double		GetLifetime()		const		{ return Lifetime;		}	/*!< Return literature lifetime */
		double		GetDnUnc()		const		{ return DnUncertainty;		}	/*!< Return negative uncertainty */
		double		GetUpUnc()		const		{ return UpUncertainty;		}	/*!< Return positive uncertainty */

	private:
		int		StateIndex;		/*!< State index */
		double		Lifetime;		/*!< Literature lifetime */
		double		DnUncertainty;		/*!< Negative uncertainty */
		double		UpUncertainty;		/*!< Positive uncertainty */

};

///
///	\class LitBranchingRatio
///
///	\brief Holder for literature branching ratio information for fitting
///
class LitBranchingRatio{		// Literature branching ratio information holder

	public:
		LitBranchingRatio()	{;}
		///
		///	Construct literature branching ratio with symmetric uncertainties:\n
		///	Initial state index = i1 \n
		///	Final state (1) index = f1 \n
		///	Final state (2) index = f2 \n
		///	Branching ratio (i1->f1)/(i1->f2) = br \n
		///	Uncertainty = e 
		///
		LitBranchingRatio(int i1, int f1, int f2, double br, double e)			{ StateIndex_I1 = i1; StateIndex_F1 = f1; StateIndex_F2 = f2; BranchingRatio = br; UpUncertainty =  e; DnUncertainty =  e; 	}
		///
		///	Construct literature branching ratio with symmetric uncertainties:\n
		///	Initial state index = i1 \n
		///	Final state (1) index = f1 \n
		///	Final state (2) index = f2 \n
		///	Branching ratio (i1->f1)/(i1->f2) = br \n
		///	Positive uncertainty = ue 
		///	Negative uncertainty = de
		///
		LitBranchingRatio(int i1, int f1, int f2, double br, double ue, double de)	{ StateIndex_I1 = i1; StateIndex_F1 = f1; StateIndex_F2 = f2; BranchingRatio = br; UpUncertainty = ue; DnUncertainty = de; 	}		// State index, lifetime, down error, up error
		~LitBranchingRatio()	{;}
		LitBranchingRatio(const LitBranchingRatio& lb);				/*!< Copy constructor */
		LitBranchingRatio& operator = (const LitBranchingRatio& lb);   		/*!< Assignment operator */

		double		GetInitialIndex() 	const		{ return StateIndex_I1;		}	/*!< Return initial state index */
		double		GetFinalIndex_1() 	const		{ return StateIndex_F1;		}	/*!< Return final state (1) */
		double		GetFinalIndex_2() 	const		{ return StateIndex_F2;		}	/*!< Return final state (2) */
		double		GetBranchingRatio() 	const		{ return BranchingRatio;	}	/*!< Return branching ratio, (i->f1)/(i->f2) */
		double		GetDnUnc() 		const		{ return DnUncertainty;		}	/*!< Return negative uncertainty */
		double		GetUpUnc() 		const		{ return UpUncertainty;		}	/*!< Return positive uncertainty */

	private:
		int		StateIndex_I1;		/*!< Initial state index */
		int		StateIndex_F1;		/*!< Final state (1) index */
		int		StateIndex_F2;		/*!< Final state (2) index */
		double		BranchingRatio;		/*!< Literature branching ratio (i->f1)/(i->f2) */
		double		DnUncertainty;		/*!< Negative uncertainty */
		double		UpUncertainty;		/*!< Positive uncertainty */

};

///
///	\class LitMixingRatio
///
///	\brief Holder for literature mixing ratio information for fitting
///
class LitMixingRatio{			// Literature mixing ratio information holder

	public:
		LitMixingRatio()	{;}
		///
		///	Construct literature mixing ratio with symmetric uncertainties:\n
		///	Initial state index = i \n
		///	Final state index = f \n
		///	Mixing ratio = d \n
		///	Uncertainty = e 
		///
		LitMixingRatio(int i, int f, double d, double e)				{ StateIndex_I = i; StateIndex_F = f; MixingRatio = d; UpUncertainty = 	e; DnUncertainty =  e;	}
		///
		///	Construct literature mixing ratio with symmetric uncertainties:\n
		///	Initial state index = i \n
		///	Final state index = f \n
		///	Mixing ratio = d \n
		///	Positive uncertainty = ue 
		///	Negative uncertainty = de
		///
		LitMixingRatio(int i, int f, double d, double ue, double de)			{ StateIndex_I = i; StateIndex_F = f; MixingRatio = d; UpUncertainty = ue; DnUncertainty = de;	}
		~LitMixingRatio()	{;}
		LitMixingRatio(const LitMixingRatio& lm);			 	/*!< Copy constructor */
		LitMixingRatio& operator = (const LitMixingRatio& lm);          	/*!< Assignment operator */

		double		GetInitialIndex()	const		{ return StateIndex_I;		}	/*!< Return initial state index */
		double		GetFinalIndex()		const		{ return StateIndex_F;		}	/*!< Return final state index */
		double		GetMixingRatio()	const		{ return MixingRatio;		}	/*!< Return mixing ratio */
		double		GetDnUnc()		const		{ return DnUncertainty;		}	/*!< Return negative uncertainty */
		double		GetUpUnc()		const		{ return UpUncertainty;		}	/*!< Return positive uncertainty */

	private:
		int 		StateIndex_I;	/*!< Initial state index */
		int 		StateIndex_F;	/*!< Final state index */
		double		MixingRatio;	/*!< Literature mixing ratio */
		double  	DnUncertainty;	/*!< Negative uncertainty */
		double		UpUncertainty;	/*!< Positive uncertainty */


};

///
///	\class LitMatrixElement
///
///	\brief Holder for literature matrix element information for fitting
///
class LitMatrixElement{

	public:
		LitMatrixElement()	{;}
		///
		///	Construct literature mixing ratio with symmetric uncertainties:\n
		///	Multipolarity = mult \n
		///	Initial state index = i \n
		///	Final state index = f \n
		///	Matrix element = me \n
		///	Uncertainty = e 
		///
		LitMatrixElement(int mult, int i, int f, double me, double e)			{ Multipolarity = mult; StateIndex_I = i; StateIndex_F = f; MatrixElement = me; UpUncertainty =	e; DnUncertainty =  e;	}
		///
		///	Construct literature mixing ratio with symmetric uncertainties:\n
		///	Multipolarity = mult \n
		///	Initial state index = i \n
		///	Final state index = f \n
		///	Matrix element = me \n
		///	Positive uncertainty = ue 
		///	Negative uncertainty = de
		///
		LitMatrixElement(int mult, int i, int f, double me, double ue, double de)	{ Multipolarity = mult; StateIndex_I = i; StateIndex_F = f; MatrixElement = me; UpUncertainty = ue; DnUncertainty = de;	}
		~LitMatrixElement()	{;}
		LitMatrixElement(const LitMatrixElement& lm);			 	/*!< Copy constructor */
		LitMatrixElement& operator = (const LitMatrixElement& lm);          	/*!< Assignment operator */

		double		GetMultipolarity()	const		{ return Multipolarity;		} 	/*!< Return multipolarity of the matrix element */
		double		GetInitialIndex()	const		{ return StateIndex_I;		}	/*!< Return initial state index */
		double		GetFinalIndex()		const		{ return StateIndex_F;		}	/*!< Return final state index */
		double		GetMatrixElement()	const		{ return MatrixElement;		}	/*!< Return mixing ratio */
		double		GetDnUnc()		const		{ return DnUncertainty;		}	/*!< Return negative uncertainty */
		double		GetUpUnc()		const		{ return UpUncertainty;		}	/*!< Return positive uncertainty */

	private:
		double		Multipolarity;	/*!< Matrix element multipolarity */
		int 		StateIndex_I;	/*!< Initial state index */
		int 		StateIndex_F;	/*!< Final state index */
		double		MatrixElement;	/*!< Literature matrix element */
		double  	DnUncertainty;	/*!< Negative uncertainty */
		double		UpUncertainty;	/*!< Positive uncertainty */
};
#endif

