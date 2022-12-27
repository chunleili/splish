#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"
#include "NonNewton.h"

using namespace SPH;
using namespace GenParam;

int NonNewton::VISCOSITY_COEFFICIENT_BOUNDARY = -1;
int NonNewton::VISCOSITY_COEFFICIENT = -1;
int NonNewton::NON_NEWTON_METHOD = -1;
int NonNewton::ENUM_NEWTONIAN = -1;
int NonNewton::ENUM_POWER_LAW = -1;
int NonNewton::POWER_INDEX = -1;
int NonNewton::CONSISTENCY_INDEX = -1;

NonNewton::NonNewton(FluidModel *model) :
NonPressureForceBase(model)
{
	std::cout<<"constructor\n";
	m_boundaryViscosity = 0.01f;
	power_index = 0.5;
	consistency_index = 1.0;

	unsigned int numParticles = m_model->numActiveParticles();

	m_strainRate.resize(numParticles, Vector6r::Zero());
	model->addField({ "strainRate", FieldType::Vector6, [&](const unsigned int i) -> Vector6r* { return &m_strainRate[i]; }, true });

	m_viscosity.resize(numParticles, 0.0);
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

void NonNewton::initParameters()
{
	VISCOSITY_COEFFICIENT_BOUNDARY = createNumericParameter("viscosityBoundary", "Viscosity coefficient (Boundary)", &m_boundaryViscosity);
	setGroup(VISCOSITY_COEFFICIENT_BOUNDARY, "Viscosity");
	setDescription(VISCOSITY_COEFFICIENT_BOUNDARY, "Coefficient for the viscosity force computation at the boundary.");

	NON_NEWTON_METHOD = createEnumParameter("nonNewtonMethod", "nonNewtonMethod", &m_nonNewtonMethod);
	setGroup(NON_NEWTON_METHOD, "Viscosity");
	setDescription(NON_NEWTON_METHOD, "Method for nonNewton.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(NON_NEWTON_METHOD));
	enumParam->addEnumValue("Newtonian", ENUM_NEWTONIAN);
	enumParam->addEnumValue("Power Law", ENUM_POWER_LAW);

	POWER_INDEX = createNumericParameter("power_index", "power_index", &power_index);
	setGroup(POWER_INDEX, "Viscosity");
	CONSISTENCY_INDEX = createNumericParameter("consistency_index", "consistency_index", &consistency_index);
	setGroup(CONSISTENCY_INDEX, "Viscosity");
}

void NonNewton::step()
{
	static int steps{0};
	std::cout<<"\nstep: "<<steps<<"\n";
	calcStrainRate();
	computeViscosity();
	std::cout<<"m_strainRate["<<100<<"]\n"<<m_strainRate[100]<<"\n";
	std::cout<<"m_strainRate["<<100<<"].norm(): "<<m_strainRate[100].norm()<<"\n";
	steps++;

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real density0 = m_model->getValue<Real>(FluidModel::DENSITY0);

	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = (static_cast<Real>(1.0) / h);

	for (int i = 0; i < (int)numParticles; i++)
	{
		const Vector3r &xi = m_model->getPosition(i);
		const Vector3r &vi = m_model->getVelocity(i);
		Vector3r &ai = m_model->getAcceleration(i);
		const Real density_i = m_model->getDensity(i);

		//////////////////////////////////////////////////////////////////////////
		// Fluid
		//////////////////////////////////////////////////////////////////////////
		forall_fluid_neighbors(
			const Vector3r &vj = fm_neighbor->getVelocity(neighborIndex);

			// Viscosity
			const Real density_j = fm_neighbor->getDensity(neighborIndex);
			ai -= invH * m_viscosity[i] * (fm_neighbor->getMass(neighborIndex) / density_j) * (vi - vj) * sim->W(xi - xj););

		//////////////////////////////////////////////////////////////////////////
		// Boundary
		//////////////////////////////////////////////////////////////////////////
		if (m_boundaryViscosity != 0.0)
		{
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
					const Vector3r a = -invH * m_boundaryViscosity * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi - vj) * sim->W(xi - xj);
					ai += a;
					bm_neighbor->addForce(xj, -m_model->getMass(i) * a););
			}
		}
	}
}


void NonNewton::reset()
{
	std::cout<<"reset\n";
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

void NonNewton::computeViscosity() 
{
	unsigned int numParticles = m_model->numActiveParticles();
	for (unsigned int i = 0; i < numParticles; ++i)
	{
		// NonNewton viscosity
		m_viscosity[i] = consistency_index * pow(m_strainRate[i].norm(), power_index - 1);
	}
}