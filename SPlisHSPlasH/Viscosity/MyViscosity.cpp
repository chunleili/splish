#include "MyViscosity.h"

namespace SPH{

MyViscosity::MyViscosity(FluidModel *model) :
	ViscosityBase(model)
{

}

MyViscosity::~MyViscosity(void)
{

}

void MyViscosity::initParameters()
{
	ViscosityBase::initParameters();


}

void MyViscosity::step()
{

}

void MyViscosity::reset()
{

}

}