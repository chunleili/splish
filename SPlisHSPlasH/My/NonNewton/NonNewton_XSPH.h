#pragma once

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "NonPressureForceBase.h"

namespace SPH
{
	/** \brief 这是修改自XSPH的粘性算法。与原有算法的区别在于它的粘性受到非牛顿的控制，并且是个vector。
     * 
	*/
	class NonNewton_XSPH : public NonPressureForceBase
	{
	protected:
        std::vector<Real> m_viscosity;
		std::vector<Real> m_boundaryViscosity;

		virtual void initParameters();

	public:
		static int VISCOSITY_COEFFICIENT_BOUNDARY;

		NonNewton_XSPH(FluidModel *model);
		virtual ~NonNewton_XSPH(void);

		virtual void step();
		virtual void reset();

		static NonPressureForceBase* creator(FluidModel* model) { return new NonNewton_XSPH(model); }
	};
}

