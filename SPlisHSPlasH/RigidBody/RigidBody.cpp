#include "RigidBody.h"
#include "SPlisHSPlasH/TimeManager.h"
using namespace SPH;

RigidBody::RigidBody(FluidModel *model):
m_model(model)
{
}

RigidBody::~RigidBody()
{
}

void RigidBody::step()
{
    computeBarycenter();
    addForce();
    translation();
    addTorque();
    rotation();
    animateParticles();
}



void RigidBody::translation()
{
    const Real h = TimeManager::getCurrent()->getTimeStepSize();
    velocity += force / total_mass * h;
    barycenter += velocity * h;
}

void RigidBody::addForce()
{
    const Vector3r g{0.0, -9.8, 0.0};
    force = g;
}

void RigidBody::animateParticles()
{
    const  int numParticles = (int) m_model->numActiveParticles();
    if (numParticles == 0)
		return;
    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static) nowait 
        for (int i = 0; i < numParticles; i++)
        { 
            Vector3r &xi = m_model->getPosition(i);
            xi += barycenter;
        }
    }
}

void RigidBody::computeBarycenter()
{
    const  int numParticles = (int) m_model->numActiveParticles();
    if (numParticles == 0)
		return;

    Vector3r sum_mass_pos = Vector3r(0.0, 0.0, 0.0);
    Real sum_mass = 0.0;

    #pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static) nowait 
        for (int i = 0; i < numParticles; i++)
        { 
            const Vector3r &xi = m_model->getPosition(i);
            const Real mass = m_model->getMass(i);
            sum_mass += mass ;
            sum_mass_pos += mass * xi;
        }
        total_mass = sum_mass;
        barycenter = sum_mass_pos / sum_mass;
    }
}