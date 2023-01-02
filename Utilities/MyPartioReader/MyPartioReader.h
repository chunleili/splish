#pragma once

#include "SPlisHSPlasH/Common.h"
#include <vector>

namespace Utilities
{
	/** \brief Class for reading and writing partio files.
	*/
	class MyPartioReader
	{
	public:
		static bool readParticles(const std::string &fileName, const Vector3r &translation, const Matrix3r &rotation, const Real scale,
			std::vector<Vector3r> &pos, std::vector<Vector3r> &vel, std::vector<Vector3r> &uv,
			std::vector<Vector3r> &normal);
	};

}
