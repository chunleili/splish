#include <math.h>
#include "RigidBody.h"
#include "SPlisHSPlasH/TimeManager.h"
#include <Eigen/Dense>

#define _USE_MATH_DEFINES

using namespace SPH;
using namespace std;

RigidBody::RigidBody(FluidModel* model):
NonPressureForceBase(model)
{
    
}

void RigidBody::setStates()
{
    unsigned int numParticles = m_model->numActiveParticles();
    for(unsigned int i=0; i< numParticles ; i++)
    {
        m_model->setParticleState(i, ParticleState::RigidBody);
    }
}

RigidBody::~RigidBody()
{
}

void RigidBody::step()
{   
    setStates();
    computeBarycenter();
    addForce();
    translation();
    // addTorque();
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

//用四元数实现刚体旋转
void RigidBody::rotation()
{   
    const  int numParticles = (int) m_model->numActiveParticles();
    if (numParticles == 0)
		return;
    TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();
    
    //转轴和转角
    Vector3r axis{0.0, 0.0, 1.0};
    axis = axis.normalized();
    float theta = 0.1  / 180.0 * 3.14159265358979323846;

    for (int i = 0; i < numParticles; i++)
    { 
        if (m_model->getParticleState(i) == ParticleState::RigidBody)
        {
            Vector3r &point = m_model->getPosition(i);
            Quaternionr p{0.0, point[0], point[1], point[2]};

            Quaternionr q{cos(theta/2), sin(theta/2) * axis[0], sin(theta/2) * axis[1], sin(theta/2) * axis[2]};
            p = q * p * q.inverse();
            point=p.vec();
        }
    }

}

void RigidBody::animateParticles()
{
    const  int numParticles = (int) m_model->numActiveParticles();
    if (numParticles == 0)
		return;

    TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

 
    for (int i = 0; i < numParticles; i++)
    { 
        if (m_model->getParticleState(i) == ParticleState::RigidBody)
        {
            Vector3r &xi = m_model->getPosition(i);
            xi += h * velocity; //rigid body translation
            // xi *= rotationMatrix;  //rigid body rotation
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