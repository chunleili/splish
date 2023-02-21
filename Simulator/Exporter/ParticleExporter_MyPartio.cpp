#include "ParticleExporter_MyPartio.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"
// #include "extern/partio/src/lib/Partio.h"
#include "extern/my_partio/Partio.h"
#include "extern/my_partio/PartioSingleton.h"

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
			writeParticlesPartio(exportFileName + ".bgeo.gz", model);
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
	auto* d = Partio::PartioSingleton::getCurrent();
	m_particleData = d->getParticlesData();

    Partio::ParticleAttribute posAttr;
    m_particleData->attributeInfo("position", posAttr);

	// print(m_particleData);
	//根据计算结果（存储在model中），更新粒子位置
	//其余的不需要更新了，想更新什么，就从model中取出来，然后更新到单例中
	std::cout << "model->numActiveParticles(): "<<model->numActiveParticles() << std::endl;
	
	// BUG FIX: 由于partio中的粒子是从bhclassic读入的，数量是固定的。当model中额外注入了粒子的时候(比如fluid block或者emitter或其他来源)，就会导致partio中的粒子数量不够用，从而导致程序崩溃。
	unsigned int numParInPartio = m_particleData->numParticles();
	unsigned int numParInModel = model->numActiveParticles();

    for (unsigned int i = 0; i < model->numActiveParticles(); i++)
    // for (unsigned int i = 0; i < m_particleData->numParticles(); i++)
    {
		int idx = i;

		//当粒子都是partio读入的时候，不需要额外添加粒子
		//当model中的粒子数量大于partio中的粒子数量时，需要额外添加粒子
		if(i >= numParInPartio)
		{
			int idx2 = m_particleData->addParticle();
			if(idx2 != idx)
				throw std::runtime_error("addParticle failed! Index mismatch!");
		}

        float* p = m_particleData->dataWrite<float>(posAttr, idx);

		const Vector3r& x = model->getPosition(i);
        p[0] = x[0];
        p[1] = x[1];
        p[2] = x[2];
    }
    Partio::write(fileName.c_str(), *m_particleData);
	std::cout << "Write the file: "<<fileName<< std::endl;
	print(m_particleData);
	std::cout << "---------------------\n\n";
}