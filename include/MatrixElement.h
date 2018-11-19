#ifndef MatrixElement_h
#define MatrixElement_h

#include <iostream>
#include <iomanip>

///
///	\class MatrixElement
///
///	\brief Holder class for matrix element information for use in fitting routines
///


class MatrixElement{

	public:
		MatrixElement()	{;}
		MatrixElement(int i, int l, int s1, int s2, double me, double mell, double meul){
									 	index 			= i; 
										lambda 			= l; 
										initialstate 		= s1; 
										finalstate 		= s2; 
										matrixElement 		= me;	
										matrixElement_ll 	= mell;
										matrixElement_ul 	= meul;
		}							/*!< Construct matrix element me, of index i, multipolarity l, between states s1 and s2 with lower and upper limits of mell and meul  */
		~MatrixElement() {;}
		MatrixElement& operator = (const MatrixElement& m);	/*!< Assignment operator */
		MatrixElement(const MatrixElement& m);			/*!< Copy constructor */

		void	SetupME(int i, int l, int s1, int s2, double me, double mell, double meul){
									 	index 			= i; 
										lambda 			= l; 
										initialstate 		= s1; 
										finalstate 		= s2; 
										matrixElement 		= me;	
										matrixElement_ll 	= mell;
										matrixElement_ul 	= meul;
		}							/*!< Define matrix element me, of index i, multipolarity l, between states s1 and s2 with lower and upper limits of mell and meul  */
		void	SetMatrixElement(double ME)				{ matrixElement = ME;		}	/*!< Define matrix element value */
		void	SetMatrixElementLowerLimit(double ME_LL)		{ matrixElement_ll = ME_LL;	}	/*!< Define matrix element lower limit */
		void	SetMatrixElementUpperLimit(double ME_UL)		{ matrixElement_ul = ME_UL;	}	/*!< Define matrix element upper limit */
		double 	GetMatrixElement() const				{ return matrixElement;		}	/*!< Return matrix element value */
		double	GetMatrixElementLowerLimit() const			{ return matrixElement_ll;	}	/*!< Return matrix element lower limit */
		double 	GetMatrixElementUpperLimit() const			{ return matrixElement_ul;	}	/*!< Return matrix element upper limit */
		int	GetIndex() const					{ return index;			}	/*!< Return matrix element index */
		int 	GetLambda() const 					{ return lambda;		}	/*!< Return matrix element mulitpolarity */
		int 	GetInitialState() const					{ return initialstate;		}	/*!< Return initial state index */
		int	GetFinalState()	const					{ return finalstate;		}	/*!< Return final state index */

		void	Print() const;		/*!< Print matrix element information */

	private:
		int 		index;			/*!< Matrix element index */
		int 		lambda;			/*!< Matrix element multipolarity */
		int 		initialstate;		/*!< Initial state (initial and final are arbitrary) */
		int 		finalstate;		/*!< Final state (initial and final are arbitrary)*/
		double		matrixElement;		/*!< Matrix element value */
		double		matrixElement_ll;	/*!< Matrix element lower limit */
		double		matrixElement_ul;	/*!< Matrix element upper limit */

};
#endif
