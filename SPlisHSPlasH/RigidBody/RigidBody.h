#ifndef __RigidBody_h__
#define __RigidBody_h__
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"
namespace SPH
{
/**
 * @brief this class is not using the PBD,
 * this class calculate the rigid body translation and rotation for particles directly
 * the merit of this way is that the particles can switch states between rigid body and SPH fluid seamlessly
 * which is useful for phase change simulation(melting and freezing)
 * 
 */
class RigidBody : NonPressureForceBase
{
private:
    Quaternionr quaternion{0.0, 0.0, 0.0, 0.0};
    Vector3r barycenter{0.0, 0.0, 0.0};
    Vector3r angular_velocity{0.0, 0.0, 0.0};
    Vector3r force{0.0, 0.0, 0.0};
    Real total_mass{0.0};

    void setStates();
    void computeBarycenter();
    void translation();
    void rotation();
    void animateParticles();
    
public:
    RigidBody(FluidModel *model);
    ~RigidBody();
    void step();
    void addForce();
    void addTorque();
    static NonPressureForceBase* creator(FluidModel* model) {return new RigidBody(model);}

    Vector3r velocity{0.0, 0.0, 0.0};
    Matrix3r rotationMatrix ;
};

}
#endif