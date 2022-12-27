#pragma once

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief NonNewton
	 * this class is based on XSPH
	 * 这个类的作用是实现非牛顿粘性。根据m_nonNewtonMethod的值来选择不同的非牛顿粘性模型。
	 * m_nonNewtonMethod = 0 时，使用牛顿粘性模型
	 * m_nonNewtonMethod = 1 时，使用幂律粘性模型
	*/
	class NonNewton : public NonPressureForceBase
	{
	private:
		int m_nonNewtonMethod = 0;
		std::vector<Vector6r> m_strainRate;
		virtual void initParameters() override;
		void calcStrainRate();
		void computeViscosityPowerLaw();
		void computeViscosityNewtonian();
	public:
		float m_boundaryViscosity;
		float m_viscosity0; //initial viscosity
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
		static int VISCOSITY0;
		NonNewton(FluidModel *model);
        virtual void init() override;
		virtual ~NonNewton(void);
        virtual void step() override final;
        virtual void reset() override final;
        static NonPressureForceBase* creator(FluidModel* model) {return new NonNewton(model);}
	};
}