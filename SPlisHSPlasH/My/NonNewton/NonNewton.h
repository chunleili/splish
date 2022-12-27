#pragma once

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
		int m_nonNewtonMethod = 0;
		std::vector<Vector6r> m_strainRate;
		virtual void initParameters() override;
		void calcStrainRate();
		void computeViscosity();
	public:
		float m_boundaryViscosity;
		std::vector<float> m_viscosity;
		float power_index;
		float consistency_index;
		static int NON_NEWTON_METHOD;
		static int VISCOSITY_COEFFICIENT;
		static int VISCOSITY_COEFFICIENT_BOUNDARY;
		static int ENUM_NEWTONIAN;
		static int ENUM_POWER_LAW;
		static int POWER_INDEX;
		static int CONSISTENCY_INDEX;
		NonNewton(FluidModel *model);
        virtual void init() override;
		virtual ~NonNewton(void);
        virtual void step() override final;
        virtual void reset() override final;
        static NonPressureForceBase* creator(FluidModel* model) {return new NonNewton(model);}
	};
}