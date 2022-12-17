#ifndef __NonNewton_h__
#define __NonNewton_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief NonNewton
	 * this class is based on XSPH
	*/
	class NonNewton : public NonPressureForceBase
	{
	private:
		Real m_viscosity;
		Real m_boundaryViscosity;
		int m_nonNewtonMethod = 0;
		std::vector<Vector6r> m_strainRate;
		virtual void initParameters() override;
		void calcStrainRate();
	public:
		static int NON_NEWTON;
		static int VISCOSITY_COEFFICIENT;
		static int VISCOSITY_COEFFICIENT_BOUNDARY;
		static int ENUM_NEWTONIAN;
		static int ENUM_SHEAR_THINNING;
		static int ENUM_SHEAR_THICKENING;
		NonNewton(FluidModel *model);
        virtual void init() override;
		virtual ~NonNewton(void);
        virtual void step() override final;
        virtual void reset() override final;
        static NonPressureForceBase* creator(FluidModel* model) {return new NonNewton(model);}
	};
}

#endif