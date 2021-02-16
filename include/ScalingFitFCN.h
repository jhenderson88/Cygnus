#ifndef ScalingFitFCN_h
#define ScalingFitFCN_h

#include <ctime>
#include <vector>

class ScalingFitFCN{ // : public ROOT::Minuit2::FCNBase{

	public:

		ScalingFitFCN()								{	ClearAll();			}
		~ScalingFitFCN()							{	;				}

		double 	up() const	 						{ return theErrorDef;			}	/*!< Required by ROOT::Minimizer */
		double 	operator()(const double*);											/*!< Required by ROOT::Minimizer */
		void	SetErrorDef(double def)						{ theErrorDef = def;			}	/*!< Required by ROOT::Minimizer */	

		void	ClearAll()
		{
			exptData.clear();
			exptUnc.clear();
			calcData.clear();
		}

		void	AddTransition(double e, double e_u, double c)
		{
			exptData.push_back(e);
			exptUnc.push_back(e_u);
			calcData.push_back(c);
		}

		void	SetData(std::vector<double> e, std::vector<double> e_u, std::vector<double> c)
		{
			exptData 	= e;
			exptUnc		= e_u;
			calcData	= c;
		}

	private:

		std::vector<double>	exptData;
		std::vector<double>	exptUnc;
		std::vector<double>	calcData;

		double				theErrorDef;

};

#endif
