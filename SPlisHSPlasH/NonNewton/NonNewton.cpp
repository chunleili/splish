#include "NonNewton.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
using namespace SPH;

int NonNewton::NON_NEWTON = -1;

NonNewton::NonNewton(FluidModel *model) :
NonPressureForceBase(model),
m_viscosity(0.0)
{
	std::cout<<"constructor\n";
	m_strainRate.resize(model->numParticles(), Vector6r::Zero());
	model->addField({ "strainRate", FieldType::Vector6, [&](const unsigned int i) -> Vector6r* { return &m_strainRate[i]; }, true });
}

NonNewton::~NonNewton(void)
{
	m_model->removeFieldByName("strainRate");
	m_strainRate.clear();
}

void NonNewton::init()
{
	std::cout<<"init\n";
	initParameters();
}

void NonNewton::step()
{
	static int steps{0};
	std::cout<<"step: "<<steps<<"\t";
	calcStrainRate();
	std::cout<<"m_strainRate["<<100<<"]"<<m_strainRate[100];
	steps++;
}

void NonNewton::reset()
{
	std::cout<<"reset\n";
	m_viscosity=0.0;
}

void NonNewton::initParameters()
{
	NON_NEWTON = createNumericParameter("nonNewtonViscosity", "nonNewton viscosity", &m_viscosity);
	setGroup(NON_NEWTON, "NonNewton");
	setDescription(NON_NEWTON, "NonNewton Params");
}

void NonNewton::calcStrainRate()
{
// shear strain rate calculation
	Simulation *sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	unsigned int numParticles = m_model->numActiveParticles();

	for (unsigned int i = 0; i < numParticles; ++i)
	{
		Vector6r strainRate;
		const Vector3r &xi = m_model->getPosition(i);
		const Vector3r &vi = m_model->getVelocity(i);
		const Real density_i = m_model->getDensity(i);

		for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
		{
			const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
			const Vector3r &xj = m_model->getPosition(neighborIndex);

			const Vector3r &vj = m_model->getVelocity(neighborIndex);

			const Vector3r gradW = sim->gradW(xi - xj);
			const Vector3r vji = vj - vi;
			const Real m = m_model->getMass(neighborIndex);
			const Real m2 = m * static_cast<Real>(2.0);
			strainRate[0] += m2 * vji[0] * gradW[0];
			strainRate[1] += m2 * vji[1] * gradW[1];
			strainRate[2] += m2 * vji[2] * gradW[2];
			strainRate[3] += m * (vji[0] * gradW[1] + vji[1] * gradW[0]);
			strainRate[4] += m * (vji[0] * gradW[2] + vji[2] * gradW[0]);
			strainRate[5] += m * (vji[1] * gradW[2] + vji[2] * gradW[1]);
		}
		strainRate = (static_cast<Real>(0.5) / density_i) * strainRate;

		m_strainRate[i] = strainRate;
		// std::cout<<"m_strainRate["<<i<<"]"<<m_strainRate[i];
		// end shear strain rate calculation
	}
}