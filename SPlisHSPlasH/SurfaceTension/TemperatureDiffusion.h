#pragma once
#ifndef __TemperatureDiffusion_h__
#define __TemperatureDiffusion_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SurfaceTensionBase.h"

namespace SPH
{
	class TemperatureDiffusion : public SurfaceTensionBase
	{
	protected:
		virtual void initParameters();
		virtual void deferredInit();
		void initValues();
		Real m_diffusivity;
		Real m_rSource;

	public:
		TemperatureDiffusion();
		TemperatureDiffusion(FluidModel* model);
		virtual ~TemperatureDiffusion(void);
		static NonPressureForceBase* creator(FluidModel* model) { return new TemperatureDiffusion(model); }
		virtual void step();
		virtual void reset();

		static int DIFFUSIVITY;
		static int R_SOURCE;
	};
}

#endif