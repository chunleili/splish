#include "MyTimeStep.h"

using namespace SPH;
using namespace std;

MyTimeStep::MyTimeStep() :
    TimeStep()
{
    cout<<"constructing MyTimeStep..."<<endl;
}

MyTimeStep::~MyTimeStep(void)
{
    cout<<"destructing MyTimeStep..."<<endl;
}

void MyTimeStep::step()
{
    cout<<"steps:"<<steps<<endl;
    steps++;
}

void MyTimeStep::resize()
{
    cout<<"resize my arrays..."<<endl;
}

void MyTimeStep::reset()
{
    steps = 0;
    cout<<"reset..."<<endl;
}