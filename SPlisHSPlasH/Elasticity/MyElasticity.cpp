#include "MyElasticity.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Utilities/Logger.h"

using namespace SPH;
using namespace GenParam;




MyElasticity::MyElasticity(FluidModel *model) :
	ElasticityBase(model)
{
}

MyElasticity::~MyElasticity(void)
{
}


void MyElasticity::initParameters()
{
    NonPressureForceBase::initParameters();
}


void MyElasticity::determineFixedParticles()
{

}

void MyElasticity::step()
{

}

void MyElasticity::reset()
{

}

