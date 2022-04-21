#include "Viscosity_XSPH.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

using namespace SPH;
using namespace GenParam;

int Viscosity_XSPH::VISCOSITY_COEFFICIENT_BOUNDARY = -1;
int Viscosity_XSPH::MU_C = -1;
int Viscosity_XSPH::TAU_C = -1;
int Viscosity_XSPH::LAMBDA = -1;


Viscosity_XSPH::Viscosity_XSPH(FluidModel *model) :
	ViscosityBase(model)
{
	m_boundaryViscosity = 0.0;

	m_muC = 0.00298;
	m_tauC = 0.02876;
	m_lambda = 4.020;

	m_strainRate.resize(model->numParticles(), Vector6r::Zero());
	model->addField({ "strain rate", FieldType::Vector6, [&](const unsigned int i) -> Real* { return &m_strainRate[i][0]; } });

	m_cassonViscosity.resize(model->numParticles(), 0.0);
	model->addField({ "strain rate", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_cassonViscosity[i]; } });
}

Viscosity_XSPH::~Viscosity_XSPH(void)
{
	m_model->removeFieldByName("strain rate");
	m_strainRate.clear();
}

void Viscosity_XSPH::initParameters()
{
	ViscosityBase::initParameters();

	VISCOSITY_COEFFICIENT_BOUNDARY = createNumericParameter("viscosityBoundary", "Viscosity coefficient (Boundary)", &m_boundaryViscosity);
	setGroup(VISCOSITY_COEFFICIENT_BOUNDARY, "Viscosity");
	setDescription(VISCOSITY_COEFFICIENT_BOUNDARY, "Coefficient for the viscosity force computation at the boundary.");
	GenParam::RealParameter* rparam = static_cast<RealParameter*>(getParameter(VISCOSITY_COEFFICIENT_BOUNDARY));
	rparam->setMinValue(0.0);

	MU_C = createNumericParameter("mu_c", "mu_c: coefficient of modified Casson model", &m_muC);
	setGroup(MU_C, "coagualtion");
	setDescription(MU_C, "mu_c: coefficient of modified Casson model");
	rparam = static_cast<GenParam::RealParameter*>(getParameter(MU_C));
	rparam->setMinValue(0.0);

	TAU_C = createNumericParameter("tau_C", "tau_C: coefficient of modified Casson model", &m_tauC);
	setGroup(TAU_C, "coagualtion");
	rparam = static_cast<GenParam::RealParameter*>(getParameter(TAU_C));
	rparam->setMinValue(0.0);

	LAMBDA = createNumericParameter("lambda", "lambda: coefficient of modified Casson model", &m_lambda);
	setGroup(LAMBDA, "coagualtion");
	setDescription(LAMBDA, "lambda: coefficient of modified Casson model");
	rparam = static_cast<GenParam::RealParameter*>(getParameter(LAMBDA));
	rparam->setMinValue(0.0);
}

void Viscosity_XSPH::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real density0 = m_model->getValue<Real>(FluidModel::DENSITY0);

	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = (static_cast<Real>(1.0) / h);

	// Compute viscosity forces (XSPH)
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &vi = m_model->getVelocity(i);
			Vector3r &ai = m_model->getAcceleration(i);
			const Real density_i = m_model->getDensity(i);


			//shear strain rate calculation
			Vector6r &strainRate = m_strainRate[i];
			
			for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
			{
				const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
				const Vector3r& xj = m_model->getPosition(neighborIndex);

				const Vector3r& vj = m_model->getVelocity(neighborIndex);

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
			//end shear strain rate calculation

			//Modified Casson Viscosity
			m_cassonViscosity[i] = sqrt(m_muC) + sqrt(m_tauC) / (sqrt(m_lambda) + sqrt(strainRate.norm()));		
			m_viscosity = m_cassonViscosity[i];
			//end: Modified Casson Viscosity

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Vector3r &vj = fm_neighbor->getVelocity(neighborIndex);

				// Viscosity
				const Real density_j = fm_neighbor->getDensity(neighborIndex);
				ai -= invH * m_viscosity * (fm_neighbor->getMass(neighborIndex) / density_j) * (vi - vj) * sim->W(xi - xj);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (m_boundaryViscosity != 0.0)
			{
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					forall_boundary_neighbors(
						const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
						const Vector3r a = -invH * m_boundaryViscosity * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi-vj)* sim->W(xi - xj);
						ai += a;
						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						Vector3r vj;
						bm_neighbor->getPointVelocity(xi, vj);
						const Vector3r a = -invH * m_boundaryViscosity * (density0 / density_i) * (vi-vj)* rho;
						ai += a;
						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				{
					forall_volume_maps(
						Vector3r vj;
						bm_neighbor->getPointVelocity(xj, vj);
						const Vector3r a = -invH * m_boundaryViscosity * (density0 * Vj / density_i) * (vi-vj)* sim->W(xi - xj);
						ai += a;
						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
					);
				}
			}
		}
	}
}


void Viscosity_XSPH::reset()
{
}

