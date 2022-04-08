#pragma once
#ifndef __MyViscosity_h__
#define __MyViscosity_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "ViscosityBase.h"

namespace SPH
{
	class MyViscosity : public ViscosityBase
	{
	protected:
		virtual void initParameters();

	public:
		MyViscosity(FluidModel* model);
		virtual ~MyViscosity(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new MyViscosity(model); }

		virtual void step();
		virtual void reset();
	};
}

#endif