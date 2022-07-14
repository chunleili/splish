#pragma once
#ifndef __Coagulation_h__
#define __Coagulation_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/SurfaceTension/SurfaceTensionBase.h"

namespace SPH
{
	class Coagulation : public SurfaceTensionBase
	{
	protected:
		virtual void initParameters();

		Real m_thresLow;
		Real m_thresHigh;
		Real m_diffusivity;
		Real m_rSource;

		Real m_pointSrcVal;
		int m_pointSrcPos=-1;

		Vector3r m_coaguBoxMin;
		Vector3r m_coaguBoxMax;

		std::vector<Real> m_ccf;

		virtual void deferredInit();

		void initValues();

	public:
		static int THRES_LOW;
		static int THRES_HIGH;
		static int DIFFUSIVITY;
		static int R_SOURCE;

		static int POINT_SRC_VAL;
		static int POINT_SRC_POS;

		static int COAGU_BOX_MIN;
		static int COAGU_BOX_MAX;

		Coagulation();
		Coagulation(FluidModel* model);
		virtual ~Coagulation(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new Coagulation(model); }

		virtual void step();
		virtual void reset();


		virtual void performNeighborhoodSearchSort();

		Real& getCcf(const unsigned int i)
		{
			return m_ccf[i];
		}

	};
}

#endif