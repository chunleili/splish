#include "ParticleExporter_MyPartio.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"
// #include "extern/partio/src/lib/Partio.h"
#include "extern/my_partio/Partio.h"

using namespace SPH;
using namespace Utilities;

ParticleExporter_MyPartio::ParticleExporter_MyPartio(SimulatorBase *base) :
	ExporterBase(base)
{
}

ParticleExporter_MyPartio::~ParticleExporter_MyPartio(void)
{
}

void ParticleExporter_MyPartio::init(const std::string& outputPath)
{
	m_exportPath = FileSystem::normalizePath(outputPath + "/mypartio");
}

void ParticleExporter_MyPartio::step(const unsigned int frame)
{
	if (!m_active)
		return;

	Simulation* sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel* model = sim->getFluidModel(i);
		std::string fileName = "ParticleData";
		if (!m_base->getValue<bool>(SimulatorBase::EXPORT_OBJECT_SPLITTING))
		{
			fileName = fileName + "_" + model->getId() + "_" + std::to_string(frame);
			std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName);
			writeParticlesPartio(exportFileName + ".bgeo", model);
		}
		else
		{
			// object splitting
			for (auto j = 0u; j < m_base->getLastObjectId(); j++)
			{
				std::string fileName2 = fileName + "_" + model->getId() + "_" + std::to_string(j) + "_" + std::to_string(frame);
				std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName2);
				writeParticlesPartio(exportFileName + ".bgeo.gz", model, j);
			}
		}
	}
}

void ParticleExporter_MyPartio::reset()
{
}

void ParticleExporter_MyPartio::setActive(const bool active)
{
	ExporterBase::setActive(active);
	if (m_active)
		FileSystem::makeDirs(m_exportPath);
}


void ParticleExporter_MyPartio::writeParticlesPartio(const std::string& fileName, FluidModel* model, const unsigned int objId)
{
    m_particleData = Partio::create();
    Partio::ParticleAttribute posAttr, densityAttr,uvAttr;
    posAttr = m_particleData->addAttribute("position", Partio::VECTOR, 3);
    densityAttr = m_particleData->addAttribute("density", Partio::FLOAT, 1);
    uvAttr = m_particleData->addAttribute("uv", Partio::VECTOR, 2);

    for (int i = 0; i < model->numActiveParticles(); i++)
    {
        int idx = m_particleData->addParticle();
        float* p = m_particleData->dataWrite<float>(posAttr, idx);
        float* den = m_particleData->dataWrite<float>(densityAttr, idx);
        float* uv = m_particleData->dataWrite<float>(uvAttr, idx);

		const Vector3r& x = model->getPosition(i);
        p[0] = x[0];
        p[1] = x[1];
        p[2] = x[2];

		den[0] = model->getDensity(i);
    }

    Partio::write(fileName.c_str(), *m_particleData);
    m_particleData->release();
}