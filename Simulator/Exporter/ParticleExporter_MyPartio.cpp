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


	/* -------------------------------------------------------------------------- */
	/*                               导出splish新增自定义场                        */
	/* -------------------------------------------------------------------------- */
	// ADD FEATURE 2023.3.1 : 增加导出splish内部有，但bhclassic中本身没有的量
	// add attributes
	std::vector<std::string> attributes;
	StringTools::tokenize(m_base->getValue<std::string>(SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES), attributes, ";");

	std::cout<<"There are "<<attributes.size()<<" user defined attributes to be exported.\n They are: \n";
	for(auto attr : attributes)
	{
		std::cout<<attr<<"\n";
	}

	std::cout<<"There are "<<model->numberOfFields()<<" fields in fluid model.\n They are: \n";
	for (unsigned int j = 0; j < model->numberOfFields(); j++)
	{	
		const FieldDescription& field = model->getField(j);
		std::cout<<field.name<<"\n";
	}

	std::cout<<"Begin to search and match for user defined attributes in fluid model.\n";
	// 搜索所有场的名字，看是否和attributes中的名字相同，相同的话，就把场的名字和场的index存入map中，然后为partio中的粒子添加属性
	std::map<unsigned int, int> attrMap;
	std::map<unsigned int, Partio::ParticleAttribute> partioAttrMap;
	for (unsigned int i = 0; i < attributes.size(); i++)
	{
		// position is exported anyway
		if (attributes[i] == "position")
		{
			attrMap[i] = -1;
			continue;
		}
		bool found = false;
		for (unsigned int j = 0; j < model->numberOfFields(); j++)
		{	
			const FieldDescription& field = model->getField(j);

			if (field.name == attributes[i])
			{
				found = true;
				std::cout<<"We found attribute "<<attributes[i]<<" in fluid model.\n It's type is: ";

				if (field.type == Scalar)
				{
					std::cout<<"Scalar\n";
					std::cout<<"Adding new attribute "<<attributes[i]<<" to partio.\n";
					attrMap[i] = j;
					partioAttrMap[i] = m_particleData->addAttribute(attributes[i].c_str(), Partio::FLOAT, 1);
				}
				else if (field.type == UInt)
				{
					std::cout<<"UInt\n";
					attrMap[i] = j;
					std::cout<<"Adding new attribute "<<attributes[i]<<" to partio.\n";
					partioAttrMap[i] = m_particleData->addAttribute(attributes[i].c_str(), Partio::INT, 1);
				}
				else if (field.type == Vector3)
				{
					std::cout<<"Vector3\n";
					std::cout<<"Adding new attribute "<<attributes[i]<<" to partio.\n";
					attrMap[i] = j;
					partioAttrMap[i] = m_particleData->addAttribute(attributes[i].c_str(), Partio::VECTOR, 3);
				}
				else
				{
					attrMap[i] = -1;
					LOG_WARN << "Only scalar and vector fields are currently supported by the partio exporter.";
				}
				break;
			}
		}
		if (!found)
		{
			attrMap[i] = -1;
			LOG_WARN << "Unknown field cannot be exported in partio file: " << attributes[i];
		}
	}

	//然后根据搜索到的属性，写出场
	for (unsigned int i = 0; i < numParInModel; i++)
	{
		for (unsigned int j = 0; j < attributes.size(); j++)
		{
			const int fieldIndex = attrMap[j];
			if (fieldIndex != -1)
			{
				const FieldDescription& field = model->getField(fieldIndex);
				if (field.type == FieldType::Scalar)
				{
					float* val = m_particleData->dataWrite<float>(partioAttrMap[j], i);
					*val = (float)*((Real*)field.getFct(i));
				}
				else if (field.type == FieldType::UInt)
				{
					int* val = m_particleData->dataWrite<int>(partioAttrMap[j], i);
					*val = (int)*((unsigned int*)field.getFct(i));
				}
				else if (field.type == FieldType::Vector3)
				{
					float* val = m_particleData->dataWrite<float>(partioAttrMap[j], i);
					Eigen::Map<Vector3r> vec((Real*)field.getFct(i));
					val[0] = (float)vec[0];
					val[1] = (float)vec[1];
					val[2] = (float)vec[2];
				}
			}
		}
	}
	/* -------------------------------------------------------------------------- */
	/*                          End of  导出splish新增自定义场                      */
	/* -------------------------------------------------------------------------- */



    Partio::write(fileName.c_str(), *m_particleData);
	std::cout << "Write partio data to the file: "<<fileName<< std::endl;
	print(m_particleData);
	std::cout << "---------------------\n\n";
}