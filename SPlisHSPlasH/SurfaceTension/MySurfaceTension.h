#pragma once
#ifndef __MySurfaceTension_h__
#define __MySurfaceTension_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SurfaceTensionBase.h"

namespace SPH
{
	class MySurfaceTension : public SurfaceTensionBase
	{
	protected:
		virtual void initParameters();

	public:
		MySurfaceTension(FluidModel* model);
		virtual ~MySurfaceTension(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new MySurfaceTension(model); }

		virtual void step();
		virtual void reset();
	};
}

#endif