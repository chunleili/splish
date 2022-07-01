#ifndef __NonNewton_h__
#define __NonNewton_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief NonNewton
	 * this class is registered in vorticity
	*/
	class NonNewton : public NonPressureForceBase
	{
	private:
		Real m_viscosity;
		virtual void initParameters() override;
	public:
		static int NON_NEWTON;
		NonNewton(FluidModel *model);
		virtual ~NonNewton(void);
        virtual void step() override;
        virtual void init() override;
        virtual void reset() override;
        static NonPressureForceBase* creator(FluidModel* model) {return new NonNewton(model);}
	};
}

#endif