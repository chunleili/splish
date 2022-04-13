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

		Real m_thresLow;
		Real m_thresHigh;
		Real m_diffusivity;
		Real m_rSource;

		std::vector<Real> m_ccf;

		virtual void deferredInit();

		void initValues();

	public:
		static int THRES_LOW;
		static int THRES_HIGH;
		static int DIFFUSIVITY;
		static int R_SOURCE;


		MySurfaceTension();
		MySurfaceTension(FluidModel* model);
		virtual ~MySurfaceTension(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new MySurfaceTension(model); }

		virtual void step();
		virtual void reset();


		virtual void performNeighborhoodSearchSort();

		Real& getCcf(const unsigned int i)
		{
			return m_ccf[i];
		}

		// FORCE_INLINE const Real getCcf(const unsigned int i) const
		// {
		// 	return m_ccf[i];
		// }

		// FORCE_INLINE Real& getCcf(const unsigned int i)
		// {
		// 	return m_ccf[i];
		// }

		// FORCE_INLINE void setCcf(const unsigned int i, const Real val)
		// {
		// 	m_ccf[i] = val;
		// }
	};
}

#endif