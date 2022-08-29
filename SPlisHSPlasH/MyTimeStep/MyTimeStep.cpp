#include "MyTimeStep.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"

using namespace SPH;
using namespace std;
using namespace GenParam;

MyTimeStep::MyTimeStep() :
    TimeStep()
{
    m_counter = 0;
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

    Simulation *sim = Simulation::getCurrent();
    const unsigned int nModels = sim->numberOfFluidModels();
    TimeManager *tm = TimeManager::getCurrent();
    const Real h = tm->getTimeStepSize();
    
    performNeighborhoodSearch();
}

void MyTimeStep::resize()
{
    cout<<"resize my arrays..."<<endl;
}

void MyTimeStep::reset()
{
    steps = 0;
    cout<<"reset..."<<endl;
    m_counter = 0;
}

void MyTimeStep::performNeighborhoodSearch()
{
	// if (Simulation::getCurrent()->zSortEnabled())
	// {
	// 	if (m_counter % 500 == 0)
	// 	{
	// 		Simulation::getCurrent()->performNeighborhoodSearchSort();
	// 		m_simulationData.performNeighborhoodSearchSort();
	// 	}
	// 	m_counter++;
	// }

	// Simulation::getCurrent()->performNeighborhoodSearch();
}