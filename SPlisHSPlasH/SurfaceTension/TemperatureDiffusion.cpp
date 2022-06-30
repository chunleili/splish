#include "TemperatureDiffusion.h"
#include "SPlisHSPlasH/Simulation.h"
#include <cstdio>
#include "SPlisHSPlasH/TimeManager.h"

using namespace SPH;
using namespace GenParam;
int TemperatureDiffusion::DIFFUSIVITY = -1;
int TemperatureDiffusion::R_SOURCE = -1;

TemperatureDiffusion::TemperatureDiffusion(FluidModel* model) :
    SurfaceTensionBase(model)
{
    model->setTemperature(1501, 100.0);
    m_diffusivity = 50.0;
    m_rSource = 0.0;
}



TemperatureDiffusion::~TemperatureDiffusion(void)
{
}

void TemperatureDiffusion::initParameters()
{
    SurfaceTensionBase::initParameters();
}

void TemperatureDiffusion::step()
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
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();

    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static)  
        for (int i = 0; i < (int)numParticles; i++)
        {
            const Vector3r &xi = m_model->getPosition(i);
            Real density_i = m_model->getDensity(i);
            
            Real &temp = model->getTemperature(i);
            Real temp_sum = 0.0;
            Real temp_old = temp;

            for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
            {
                const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
                const Vector3r &xj = model->getPosition(neighborIndex);
                const Real grawWNorm = sim->gradW(xi - xj).norm();

                Real density_j = m_model->getDensity(neighborIndex);

                Real temp_j = model->getTemperature(neighborIndex);
                temp_sum += (m_diffusivity * m_model->getMass(neighborIndex) 
                / (density_j * density_i) * (temp_j - temp) * grawWNorm + m_rSource) * dt;
            }
            temp = temp_sum + temp_old;
        }
    }
}

void TemperatureDiffusion::reset()
{}

void TemperatureDiffusion::deferredInit()
{}

void TemperatureDiffusion::initValues()
{}
