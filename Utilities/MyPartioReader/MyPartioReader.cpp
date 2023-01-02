#include "MyPartioReader.h"
#include "extern/my_partio/Partio.h"
#include "../FileSystem.h"

using namespace Utilities;


bool MyPartioReader::readParticles(const std::string &fileName, const Vector3r &translation, const Matrix3r &rotation, const Real scale,
	std::vector<Vector3r> &positions, std::vector<Vector3r> &velocities, std::vector<Vector3r> &uv)
{
	if (!FileSystem::fileExists(fileName))
		return false;

	Partio::ParticlesDataMutable* data = Partio::read(fileName.c_str());

    std::cout<<"Reading Partio using MyPartioReader\n";
    print(data);
    std::cout<<"End reading Partio using MyPartioReader\n";

	if (!data)
		return false;

	unsigned int posIndex = 0xffffffff;
	unsigned int velIndex = 0xffffffff;
	unsigned int uvIndex = 0xffffffff;

	for (int i = 0; i < data->numAttributes(); i++)
	{
		Partio::ParticleAttribute attr;
		data->attributeInfo(i, attr);
		if (attr.name == "position")
			posIndex = i;
		else if (attr.name == "velocity")
			velIndex = i;
        else if (attr.name == "uv")
			uvIndex = i;
	}

	Partio::ParticleAttribute attr;

    std::cout<<"Reading uv: "<<uvIndex<<"\n!!!!!!!!!!!\n";
    std::cout<<"Reading position: "<<posIndex<<"\n!!!!!!!!!!!\n";


	if (posIndex != 0xffffffff)
	{
		unsigned int fSize = (unsigned int) positions.size();
		positions.resize(fSize + data->numParticles());
		data->attributeInfo(posIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *pos = data->data<float>(attr, i);
			Vector3r x(pos[0], pos[1], pos[2]);
			x = rotation * (x*scale) + translation;
			positions[i + fSize] = x;
		}
	}

	if (velIndex != 0xffffffff)
	{
		unsigned int fSize = (unsigned int) velocities.size();
		velocities.resize(fSize + data->numParticles());
		data->attributeInfo(velIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *vel = data->data<float>(attr, i);
			Vector3r v(vel[0], vel[1], vel[2]);
			velocities[i + fSize] = v;
		}
	}

	else
	{
		unsigned int fSize = (unsigned int) velocities.size();
		velocities.resize(fSize + data->numParticles());
		for (int i = 0; i < data->numParticles(); i++)
			velocities[i + fSize].setZero();
	}

	if (uvIndex != 0xffffffff)

	{
        std::cout<<"I am Reading uv: ";
        std::cout<<uvIndex<<"\n";

		unsigned int fSize = (unsigned int) uv.size();
		uv.resize(fSize + data->numParticles());
		data->attributeInfo(uvIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *uv_ = data->data<float>(attr, i);
			Vector3r uv__(uv_[0], uv_[1], uv_[2]);
			uv[i + fSize] = uv__;

		}
	}
    

	data->release();
	return true;
}
