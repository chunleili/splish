#pragma once
#ifndef __MyVorticity_h__
#define __MyVorticity_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/Vorticity/VorticityBase.h"

namespace SPH
{
	class MyVorticity : public VorticityBase
	{
	protected:
		Real m_vorticityCoeff;

		virtual void initParameters();

	public:
		static int VORTICITY_COEFFICIENT;

		MyVorticity(FluidModel *model);
		virtual ~MyVorticity(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new MyVorticity(model); }

		virtual void step();
		virtual void reset();
		
		virtual void performNeighborhoodSearchSort();

	};
}

#endif