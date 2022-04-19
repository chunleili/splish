#include "MySurfaceTension.h"
#include "SPlisHSPlasH\Simulation.h"
#include <cstdio>

using namespace SPH;


int MySurfaceTension::THRES_HIGH  = -1;
int MySurfaceTension::THRES_LOW = -1;
int MySurfaceTension::DIFFUSIVITY = -1;
int MySurfaceTension::R_SOURCE = -1;



MySurfaceTension::MySurfaceTension(FluidModel* model) :
    SurfaceTensionBase(model), 
    m_ccf(),
    m_thresHigh(static_cast<Real>(1.0)),
    m_thresLow(static_cast<Real>(0.1)),
    m_diffusivity(static_cast<Real>(0.1)),
    m_rSource(static_cast<Real>(0.0))
{
    m_ccf.resize(model->numParticles(), 0.0);

    model->addField({ "ccf field", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_ccf[i]; } });

}



MySurfaceTension::~MySurfaceTension(void)
{
    m_model->removeFieldByName("ccf field");
    m_ccf.clear();
}

void MySurfaceTension::initParameters()
{
    SurfaceTensionBase::initParameters();

    THRES_HIGH = createNumericParameter("m_thresHigh", "threshold high", &m_thresHigh);
	setGroup(THRES_HIGH, "coagualtion");
	setDescription(THRES_HIGH, "higher treshold for coagualtion");
	GenParam::RealParameter* rparam = static_cast<GenParam::RealParameter*>(getParameter(THRES_HIGH));
	rparam->setMinValue(0.0);

    THRES_LOW = createNumericParameter("m_thresLow", "threshold low", &m_thresLow);
    setGroup(THRES_LOW, "coagualtion");
    setDescription(THRES_LOW, "lower treshold for coagualtion");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(THRES_LOW));
    rparam->setMinValue(0.0);

    DIFFUSIVITY = createNumericParameter("m_diffusivity", "diffusivity", &m_diffusivity);
    setGroup(DIFFUSIVITY, "coagualtion");
    setDescription(DIFFUSIVITY, "diffusivity of the convection-diffusion equation");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(DIFFUSIVITY));
    rparam->setMinValue(0.0);

    R_SOURCE = createNumericParameter("m_rSource", "rSource", &m_rSource);
    setGroup(R_SOURCE, "coagualtion");
    setDescription(R_SOURCE, "source term of the convection-diffusion equation");
    rparam = static_cast<GenParam::RealParameter*>(getParameter(R_SOURCE));
    rparam->setMinValue(0.0);
}

void MySurfaceTension::step()
{
    Simulation* sim = Simulation::getCurrent();
    const unsigned int numParticles = m_model->numActiveParticles();
    const Real k = m_surfaceTension;
    const Real kb = m_surfaceTensionBoundary;
    const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
    const Real diameter = static_cast<Real>(2.0) * radius;
    const Real diameter2 = diameter * diameter;
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    const Real density0 = m_model->getDensity0();
    FluidModel* model = m_model;
    const Real h = sim->getSupportRadius();

    // Compute forces
//#pragma omp parallel default(shared)
    {
//#pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            const Vector3r &xi = m_model->getPosition(i);
            Real &ccf_i = getCcf(i);

            Real density_i = m_model->getDensity(i);

            Real source = 0.0;

            if (xi[0] >= -0.35 && xi[0] <= 0.35 && xi[1] >= -0.35 && xi[1] <= 0.35 && xi[2] >= -0.35 && xi[2] <= 0.35)
            {
                source = m_rSource;
            }
            //////////////////////////////////////////////////////////////////////////
            // Fluid
            //////////////////////////////////////////////////////////////////////////
            for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
            {
                const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
                const Vector3r &xj = model->getPosition(neighborIndex);

                Real density_j = m_model->getDensity(neighborIndex);
                Real ccf_j = getCcf(neighborIndex);
                ccf_i += m_diffusivity * m_model->getMass(neighborIndex) / (density_j * density_i) * (ccf_i - ccf_j) * sim->gradW(xi - xj).norm() + source;
            }
            ccf_i /= numParticles;
        }
    }
    determteeFixedParticles();

}

void MySurfaceTension::reset()
{

}


void MySurfaceTension::performNeighborhoodSearchSort()
{

    const unsigned int numPart = m_model->numActiveParticles();
    if (numPart == 0)
        return;

    Simulation* sim = Simulation::getCurrent();
    auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
    d.sort_field(&m_ccf[0]);

}

void MySurfaceTension::deferredInit()
{
    initValues();
}

void MySurfaceTension::initValues()
{

}


void MySurfaceTension::determteeFixedParticles()
{
    const unsigned int numParticles = m_model->numActiveParticles();

    for (int i = 0; i < (int)numParticles; i++)
        {
            const Vector3r& x = m_model->getPosition(i);
            if (m_ccf[i] > m_thresHigh)
            {
                m_model->setParticleState(i, ParticleState::Fixed);
                // printf("state is %d\n", m_model->getParticleState(i));
            }
        }
}