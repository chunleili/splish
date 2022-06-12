#include "SimulatorBase.h"
#include "Exporter/ExporterBase.h"
#include "Exporter/ParticleExporter_Partio.h"
#include "Exporter/ParticleExporter_VTK.h"
#include "Exporter/RigidBodyExporter_BIN.h"
#include "Exporter/RigidBodyExporter_OBJ.h"
#include "Exporter/RigidBodyExporter_VTK.h"
#include "Exporter/ParticleExporter_xyz.h"

using namespace SPH;

void SimulatorBase::createExporters()
{
	addParticleExporter("enablePartioExport", "Partio Exporter", "Enable/disable partio export.", new ParticleExporter_Partio(this));
	addParticleExporter("enableVTKExport", "VTK Exporter", "Enable/disable VTK export.", new ParticleExporter_VTK(this));
	addParticleExporter("enableXyzExport", "xyz Exporter", "Enable/disable xyz export.", new ParticleExporter_xyz(this));


	addRigidBodyExporter("enableRigidBodyExport", "Rigid Body Exporter", "Enable/disable rigid body BIN export.", new RigidBodyExporter_BIN(this));
	addRigidBodyExporter("enableRigidBodyOBJExport", "Rigid Body OBJ Exporter", "Enable/disable rigid body OBJ export.", new RigidBodyExporter_OBJ(this));
	addRigidBodyExporter("enableRigidBodyVTKExport", "Rigid Body VTK Exporter", "Enable/disable rigid body VTK export.", new RigidBodyExporter_VTK(this));

}