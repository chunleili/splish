#include "Viscosity_Casson.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"

using namespace SPH;
using namespace GenParam;

int Viscosity_Casson::ITERATIONS = -1;
int Viscosity_Casson::MAX_ITERATIONS = -1;
int Viscosity_Casson::MAX_ERROR = -1;
int Viscosity_Casson::VISCOSITY_COEFFICIENT_BOUNDARY = -1;
int Viscosity_Casson::VISCOSITY_FIELD = -1;
Vector3r Viscosity_Casson::m_boxMin(0., 0., 0.);
Vector3r Viscosity_Casson::m_boxMax(1., 1., 1.);

Viscosity_Casson::Viscosity_Casson(FluidModel *model) :
	ViscosityBase(model), m_vDiff()
{
	m_maxIter = 100;
	m_maxError = static_cast<Real>(0.01);
	m_iterations = 0;
	m_boundaryViscosity = 0.0;
	m_tangentialDistanceFactor = static_cast<Real>(0.5);

	m_vDiff.resize(model->numParticles(), Vector3r::Zero());

	model->addField({ "velocity difference", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vDiff[i][0]; }, true });

	m_cassonViscosity.resize(model->numParticles(), 10000.0);
	model->addField({ "casson Viscosity", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_cassonViscosity[i]; }, true });
}

Viscosity_Casson::~Viscosity_Casson(void)
{
	m_model->removeFieldByName("velocity difference");

	m_vDiff.clear();

	m_model->removeFieldByName("casson Viscosity");
	m_cassonViscosity.clear();
}

void Viscosity_Casson::initParameters()
{
	ViscosityBase::initParameters();

	VISCOSITY_COEFFICIENT_BOUNDARY = createNumericParameter("viscosityBoundary", "Viscosity coefficient (Boundary)", &m_boundaryViscosity);
	setGroup(VISCOSITY_COEFFICIENT_BOUNDARY, "Viscosity");
	setDescription(VISCOSITY_COEFFICIENT_BOUNDARY, "Coefficient for the viscosity force computation at the boundary.");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(VISCOSITY_COEFFICIENT_BOUNDARY));
	rparam->setMinValue(0.0);

	ITERATIONS = createNumericParameter("viscoIterations", "Iterations", &m_iterations);
	setGroup(ITERATIONS, "Viscosity");
	setDescription(ITERATIONS, "Iterations required by the viscosity solver.");
	getParameter(ITERATIONS)->setReadOnly(true);

	MAX_ITERATIONS = createNumericParameter("viscoMaxIter", "Max. iterations (visco)", &m_maxIter);
	setGroup(MAX_ITERATIONS, "Viscosity");
	setDescription(MAX_ITERATIONS, "Max. iterations of the viscosity solver.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS))->setMinValue(1);

	MAX_ERROR = createNumericParameter("viscoMaxError", "Max. visco error", &m_maxError);
	setGroup(MAX_ERROR, "Viscosity");
	setDescription(MAX_ERROR, "Max. error of the viscosity solver.");
	rparam = static_cast<RealParameter*>(getParameter(MAX_ERROR));
	rparam->setMinValue(static_cast<Real>(1e-6));
}



void Viscosity_Casson::matrixVecProd(const Real* vec, Real *result, void *userData)
{
	Simulation *sim = Simulation::getCurrent();
	Viscosity_Casson *visco = (Viscosity_Casson*)userData;
	FluidModel *model = visco->getModel();
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	const Real h = sim->getSupportRadius();
	const Real h2 = h*h;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const Real density0 = model->getDensity0();
	const Real mu = visco->m_viscosity * density0;
	const Real mub = visco->m_boundaryViscosity * density0;
	const Real sphereVolume = static_cast<Real>(4.0 / 3.0 * M_PI) * h2*h;

	Real d = 10.0;
	if (sim->is2DSimulation())
		d = 8.0;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = model->getPosition(i);
			Vector3r ai;
			ai.setZero();
			const Real density_i = model->getDensity(i);
			const Vector3r &vi = Eigen::Map<const Vector3r>(&vec[3 * i]);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++) 
    		{ 
        		const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j); 
        		const Vector3r &xj = model->getPosition(neighborIndex); 
				const Real density_j = model->getDensity(neighborIndex);
				const Vector3r gradW = sim->gradW(xi - xj);

				const Vector3r &vj = Eigen::Map<const Vector3r>(&vec[3 * neighborIndex]);
				const Vector3r xixj = xi - xj;

				// ai += d * mu * (model->getMass(neighborIndex) / density_j) * (vi - vj).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
				// ai += d * (mu + visco->m_cassonViscosity[neighborIndex]) * (model->getMass(neighborIndex) / density_j) * (vi - vj).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
				// if ((xj[0] > m_boxMin[0]) && (xj[1] > m_boxMin[1]) && (xj[2] > m_boxMin[2]) &&
				// 	(xj[0] < m_boxMax[0]) && (xj[1] < m_boxMax[1]) && (xj[2] < m_boxMax[2]))
				if (model->getTemperature(neighborIndex) < 1.0)
				{
					ai += d * mu * visco->m_cassonViscosity[neighborIndex] * (model->getMass(neighborIndex) / density_j) * (vi - vj).dot(xixj) / (xixj.squaredNorm() + 0.01 * h2) * gradW;
				}
				else
				{
					ai += d * mu  * (model->getMass(neighborIndex) / density_j) * (vi - vj).dot(xixj) / (xixj.squaredNorm() + 0.01 * h2) * gradW;
				}
			};
			
			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (mub != 0.0)
			{
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					forall_boundary_neighbors(
						const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
						const Vector3r xixj = xi - xj;
						const Vector3r gradW = sim->gradW(xixj);
						// const Vector3r a = d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
						Vector3r a;
						if ((xj[0] > m_boxMin[0]) && (xj[1] > m_boxMin[1]) && (xj[2] > m_boxMin[2]) &&
							(xj[0] < m_boxMax[0]) && (xj[1] < m_boxMax[1]) && (xj[2] < m_boxMax[2])) {
							a = d * mub * visco->m_cassonViscosity[neighborIndex] * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
						} else {
							a = d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
						}

						ai += a;);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						const Vector3r xixj = xi - xj;
						Vector3r normal = -xixj;
						const Real nl = normal.norm();
						if (nl > static_cast<Real>(0.0001))
						{
							normal /= nl;

							Vector3r t1;
							Vector3r t2;
							MathFunctions::getOrthogonalVectors(normal, t1, t2);

							const Real dist = visco->m_tangentialDistanceFactor * h;
							const Vector3r x1 = xj - t1 * dist;
							const Vector3r x2 = xj + t1 * dist;
							const Vector3r x3 = xj - t2 * dist;
							const Vector3r x4 = xj + t2 * dist;

							const Vector3r xix1 = xi - x1;
							const Vector3r xix2 = xi - x2;
							const Vector3r xix3 = xi - x3;
							const Vector3r xix4 = xi - x4;

							const Vector3r gradW1 = sim->gradW(xix1);
							const Vector3r gradW2 = sim->gradW(xix2);
							const Vector3r gradW3 = sim->gradW(xix3);
							const Vector3r gradW4 = sim->gradW(xix4);

							// each sample point represents the quarter of the volume inside of the boundary
							const Real vol = static_cast<Real>(0.25) * rho * sphereVolume;

 							Vector3r v1;
 							Vector3r v2;
							Vector3r v3;
							Vector3r v4;
 							bm_neighbor->getPointVelocity(x1, v1);
 							bm_neighbor->getPointVelocity(x2, v2);
							bm_neighbor->getPointVelocity(x3, v3);
							bm_neighbor->getPointVelocity(x4, v4);

 							// compute forces for both sample point
 							const Vector3r a1 = d * mub * vol * (vi).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
 							const Vector3r a2 = d * mub * vol * (vi).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
							const Vector3r a3 = d * mub * vol * (vi).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
							const Vector3r a4 = d * mub * vol * (vi).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
 							ai += a1 + a2 + a3 + a4;
 						}
					);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				{
					forall_volume_maps(
						const Vector3r xixj = xi - xj;
						Vector3r normal = -xixj;
						const Real nl = normal.norm();
						if (nl > static_cast<Real>(0.0001))
						{
							normal /= nl;

							Vector3r t1;
							Vector3r t2;
							MathFunctions::getOrthogonalVectors(normal, t1, t2);

							const Real dist = visco->m_tangentialDistanceFactor * h;
							const Vector3r x1 = xj - t1*dist;
							const Vector3r x2 = xj + t1*dist;
							const Vector3r x3 = xj - t2*dist;
							const Vector3r x4 = xj + t2*dist;

							const Vector3r xix1 = xi - x1;
							const Vector3r xix2 = xi - x2;
							const Vector3r xix3 = xi - x3;
							const Vector3r xix4 = xi - x4;

							const Vector3r gradW1 = sim->gradW(xix1);
							const Vector3r gradW2 = sim->gradW(xix2);
							const Vector3r gradW3 = sim->gradW(xix3);
							const Vector3r gradW4 = sim->gradW(xix4);

							// each sample point represents the quarter of the volume inside of the boundary
							const Real vol = static_cast<Real>(0.25) * Vj;

 							Vector3r v1;
 							Vector3r v2;
							Vector3r v3;
							Vector3r v4;
 							bm_neighbor->getPointVelocity(x1, v1);
 							bm_neighbor->getPointVelocity(x2, v2);
							bm_neighbor->getPointVelocity(x3, v3);
							bm_neighbor->getPointVelocity(x4, v4);

 							// compute forces for both sample point
 							const Vector3r a1 = d * mub * vol * (vi).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
 							const Vector3r a2 = d * mub * vol * (vi).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
							const Vector3r a3 = d * mub * vol * (vi).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
							const Vector3r a4 = d * mub * vol * (vi).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
 							ai += a1 + a2 + a3 + a4;
						}
					);
				}
			}

			result[3 * i] = vec[3 * i] - dt / density_i*ai[0];
			result[3 * i + 1] = vec[3 * i + 1] - dt / density_i*ai[1];
			result[3 * i + 2] = vec[3 * i + 2] - dt / density_i*ai[2];
		}
	}
}


#ifdef USE_BLOCKDIAGONAL_PRECONDITIONER


void Viscosity_Casson::diagonalMatrixElement(const unsigned int i, Matrix3r &result, void *userData)
{
	// Diagonal element
	Simulation *sim = Simulation::getCurrent();
	Viscosity_Casson *visco = (Viscosity_Casson*)userData;
	FluidModel *model = visco->getModel();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real density0 = model->getDensity0();

	Real d = 10.0;
	if (sim->is2DSimulation())
		d = 8.0;

	const Real h = sim->getSupportRadius();
	const Real h2 = h*h;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const Real mu = visco->m_viscosity * density0;
	const Real mub = visco->m_boundaryViscosity * density0;
	const Real sphereVolume = static_cast<Real>(4.0 / 3.0 * M_PI) * h2*h;

	const Real density_i = model->getDensity(i);

    result.setZero();
	
	const Vector3r &xi = model->getPosition(i);

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_in_same_phase(
		const Real density_j = model->getDensity(neighborIndex);
		const Vector3r gradW = sim->gradW(xi - xj);
		const Vector3r xixj = xi - xj;
		result += d * mu * (model->getMass(neighborIndex) / density_j) / (xixj.squaredNorm() + 0.01*h2) * (gradW * xixj.transpose());
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (mub != 0.0)
	{
		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		{
			forall_boundary_neighbors(
				const Vector3r xixj = xi - xj;
				const Vector3r gradW = sim->gradW(xixj);
				result += d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) / (xixj.squaredNorm() + 0.01*h2) * (gradW * xixj.transpose());
			);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		{
			forall_density_maps(
				const Vector3r xixj = xi - xj;
				Vector3r normal = -xixj;
				const Real nl = normal.norm();
				if (nl > static_cast<Real>(0.0001))
				{
					normal /= nl;

					Vector3r t1;
					Vector3r t2;
					MathFunctions::getOrthogonalVectors(normal, t1, t2);

					const Real dist = visco->m_tangentialDistanceFactor * h;
					const Vector3r x1 = xj - t1*dist;
					const Vector3r x2 = xj + t1*dist;
					const Vector3r x3 = xj - t2*dist;
					const Vector3r x4 = xj + t2*dist;

					const Vector3r xix1 = xi - x1;
					const Vector3r xix2 = xi - x2;
					const Vector3r xix3 = xi - x3;
					const Vector3r xix4 = xi - x4;

					const Vector3r gradW1 = sim->gradW(xix1);
					const Vector3r gradW2 = sim->gradW(xix2);
					const Vector3r gradW3 = sim->gradW(xix3);
					const Vector3r gradW4 = sim->gradW(xix4);

					// each sample point represents the quarter of the volume inside of the boundary
					const Real vol = static_cast<Real>(0.25) * rho * sphereVolume;

					// compute forces for both sample point
					result += d * mub * vol / (xix1.squaredNorm() + 0.01*h2) * (gradW1 * xix1.transpose());
					result += d * mub * vol / (xix2.squaredNorm() + 0.01*h2) * (gradW2 * xix2.transpose());
					result += d * mub * vol / (xix3.squaredNorm() + 0.01*h2) * (gradW3 * xix3.transpose());
					result += d * mub * vol / (xix4.squaredNorm() + 0.01*h2) * (gradW4 * xix4.transpose());
				}
			);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		{
			forall_volume_maps(
				const Vector3r xixj = xi - xj;
				Vector3r normal = -xixj;
				const Real nl = normal.norm();
				if (nl > static_cast<Real>(0.0001))
				{
					normal /= nl;

					Vector3r t1;
					Vector3r t2;
					MathFunctions::getOrthogonalVectors(normal, t1, t2);

					const Real dist = visco->m_tangentialDistanceFactor * h;
					const Vector3r x1 = xj - t1*dist;
					const Vector3r x2 = xj + t1*dist;
					const Vector3r x3 = xj - t2*dist;
					const Vector3r x4 = xj + t2*dist;

					const Vector3r xix1 = xi - x1;
					const Vector3r xix2 = xi - x2;
					const Vector3r xix3 = xi - x3;
					const Vector3r xix4 = xi - x4;

					const Vector3r gradW1 = sim->gradW(xix1);
					const Vector3r gradW2 = sim->gradW(xix2);
					const Vector3r gradW3 = sim->gradW(xix3);
					const Vector3r gradW4 = sim->gradW(xix4);

					// each sample point represents the quarter of the volume inside of the boundary
					const Real vol = static_cast<Real>(0.25) * Vj;

					// compute forces for both sample point
					result += d * mub * vol / (xix1.squaredNorm() + 0.01*h2) * (gradW1 * xix1.transpose());
					result += d * mub * vol / (xix2.squaredNorm() + 0.01*h2) * (gradW2 * xix2.transpose());
					result += d * mub * vol / (xix3.squaredNorm() + 0.01*h2) * (gradW3 * xix3.transpose());
					result += d * mub * vol / (xix4.squaredNorm() + 0.01*h2) * (gradW4 * xix4.transpose());
				}
			);
		}
	}
	result = Matrix3r::Identity() - (dt / density_i) * result;
}


#else

void Viscosity_Casson::diagonalMatrixElement(const unsigned int i, Vector3r &result, void *userData)
{
	// Diagonal element
	Simulation *sim = Simulation::getCurrent();
	Viscosity_Casson *visco = (Viscosity_Casson*)userData;
	FluidModel *model = visco->getModel();

	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	
	const Real h = sim->getSupportRadius();
	const Real h2 = h*h;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const Real density0 = model->getDensity0();
	const Real mu = visco->m_viscosity * density0;
	const Real mub = visco->m_boundaryViscosity * density0;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	Real d = 10.0;
	if (sim->is2DSimulation())
		d = 8.0;
	const Real sphereVolume = 4.0 / 3.0 * M_PI * h2*h;

	const Real density_i = model->getDensity(i);

	result.setZero();

	const Vector3r &xi = model->getPosition(i);

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_in_same_phase(
		const Real density_j = model->getDensity(neighborIndex);
		const Vector3r gradW = sim->gradW(xi - xj);
		const Vector3r xixj = xi - xj;
		Matrix3r r = d * mu * (model->getMass(neighborIndex) / density_j) / (xixj.squaredNorm() + 0.01*h2) * (gradW * xixj.transpose());
		result += r.diagonal();
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (mub != 0.0)
	{
		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		{
			forall_boundary_neighbors(
				const Vector3r xixj = xi - xj;
				const Vector3r gradW = sim->gradW(xixj);
				Matrix3r r = d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) / (xixj.squaredNorm() + 0.01*h2) * (gradW * xixj.transpose());
				result += r.diagonal();
			);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		{
			forall_density_maps(
				const Vector3r xixj = xi - xj;
				Vector3r normal = -xixj;
				const Real nl = normal.norm();
				if (nl > 0.0001)
				{
					normal /= nl;

					Vector3r t1;
					Vector3r t2;
					MathFunctions::getOrthogonalVectors(normal, t1, t2);

					const Real dist = visco->m_tangentialDistanceFactor * h;
					const Vector3r x1 = xj - t1*dist;
					const Vector3r x2 = xj + t1*dist;
					const Vector3r x3 = xj - t2*dist;
					const Vector3r x4 = xj + t2*dist;

					const Vector3r xix1 = xi - x1;
					const Vector3r xix2 = xi - x2;
					const Vector3r xix3 = xi - x3;
					const Vector3r xix4 = xi - x4;

					const Vector3r gradW1 = sim->gradW(xix1);
					const Vector3r gradW2 = sim->gradW(xix2);
					const Vector3r gradW3 = sim->gradW(xix3);
					const Vector3r gradW4 = sim->gradW(xix4);

					// each sample point represents the quarter of the volume inside of the boundary
					const Real vol = 0.25 * rho * sphereVolume;

					// compute forces for both sample point
					Matrix3r r = d * mub * vol / (xix1.squaredNorm() + 0.01*h2) * (gradW1 * xix1.transpose());
					r += d * mub * vol / (xix2.squaredNorm() + 0.01*h2) * (gradW2 * xix2.transpose());
					r += d * mub * vol / (xix3.squaredNorm() + 0.01*h2) * (gradW3 * xix3.transpose());
					r += d * mub * vol / (xix4.squaredNorm() + 0.01*h2) * (gradW4 * xix4.transpose());
					result += r.diagonal();
				}
			);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		{
			forall_volume_maps(
				const Vector3r xixj = xi - xj;
				Vector3r normal = -xixj;
				const Real nl = normal.norm();
				if (nl > 0.0001)
				{
					normal /= nl;

					Vector3r t1;
					Vector3r t2;
					MathFunctions::getOrthogonalVectors(normal, t1, t2);

					const Real dist = visco->m_tangentialDistanceFactor * h;
					const Vector3r x1 = xj - t1*dist;
					const Vector3r x2 = xj + t1*dist;
					const Vector3r x3 = xj - t2*dist;
					const Vector3r x4 = xj + t2*dist;

					const Vector3r xix1 = xi - x1;
					const Vector3r xix2 = xi - x2;
					const Vector3r xix3 = xi - x3;
					const Vector3r xix4 = xi - x4;

					const Vector3r gradW1 = sim->gradW(xix1);
					const Vector3r gradW2 = sim->gradW(xix2);
					const Vector3r gradW3 = sim->gradW(xix3);
					const Vector3r gradW4 = sim->gradW(xix4);

					// each sample point represents the quarter of the volume inside of the boundary
					const Real vol = 0.25 * Vj;

					// compute forces for both sample point
					Matrix3r r = d * mub * vol / (xix1.squaredNorm() + 0.01*h2) * (gradW1 * xix1.transpose());
					r += d * mub * vol / (xix2.squaredNorm() + 0.01*h2) * (gradW2 * xix2.transpose());
					r += d * mub * vol / (xix3.squaredNorm() + 0.01*h2) * (gradW3 * xix3.transpose());
					r += d * mub * vol / (xix4.squaredNorm() + 0.01*h2) * (gradW4 * xix4.transpose());
					result += r.diagonal();
				}
			);
		}
	}
	result = Vector3r::Ones() - (dt / density_i) * result;
}

#endif

void Viscosity_Casson::computeCassonViscosity()
{
	const int numParticles = (int) m_model->numActiveParticles();

	const Real density0 = m_model->getDensity0();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	//////////////////////////////////////////////////////////////////////////
	// locally dynamically changing viscosity
	//////////////////////////////////////////////////////////////////////////
	const Real curTime = TimeManager::getCurrent()->getTime();
	const Real startTime = 0.0;
	const Real endTime = 5.0;
	const Real endViscosity = 5000.0;

	const Real m_muC = 0.00298;
	const Real m_tauC = 0.02876;
	const Real m_lambda = 4.020;
	// shear strain rate calculation
	Simulation *sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	for (int i = 0; i < numParticles; ++i)
	{
		Vector6r strainRate;
		const Vector3r &xi = m_model->getPosition(i);
		const Vector3r &vi = m_model->getVelocity(i);
		const Real density_i = m_model->getDensity(i);

		for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i); j++)
		{
			const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);
			const Vector3r &xj = m_model->getPosition(neighborIndex);

			const Vector3r &vj = m_model->getVelocity(neighborIndex);

			const Vector3r gradW = sim->gradW(xi - xj);
			const Vector3r vji = vj - vi;
			const Real m = m_model->getMass(neighborIndex);
			const Real m2 = m * static_cast<Real>(2.0);
			strainRate[0] += m2 * vji[0] * gradW[0];
			strainRate[1] += m2 * vji[1] * gradW[1];
			strainRate[2] += m2 * vji[2] * gradW[2];
			strainRate[3] += m * (vji[0] * gradW[1] + vji[1] * gradW[0]);
			strainRate[4] += m * (vji[0] * gradW[2] + vji[2] * gradW[0]);
			strainRate[5] += m * (vji[1] * gradW[2] + vji[2] * gradW[1]);
		}
		strainRate = (static_cast<Real>(0.5) / density_i) * strainRate;
		// end shear strain rate calculation
		

		m_cassonViscosity[i] = sqrt(m_muC) + sqrt(m_tauC) / (sqrt(m_lambda) + sqrt(strainRate.norm()));
		// end: Modified Casson Viscosity


		m_cassonViscosity[i] += (curTime - startTime) / (endTime - startTime) * endViscosity;
		if (curTime > endTime)
			m_cassonViscosity[i] = endViscosity;
	}
}

void Viscosity_Casson::step()
{
	const int numParticles = (int) m_model->numActiveParticles();
	// prevent solver from running with a zero-length vector
	if (numParticles == 0)
		return;
	const Real density0 = m_model->getDensity0();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	// computeCassonViscosity();
	//////////////////////////////////////////////////////////////////////////
	// Init linear system solver and preconditioner
	//////////////////////////////////////////////////////////////////////////
	MatrixReplacement A(3*m_model->numActiveParticles(), matrixVecProd, (void*) this);
	m_solver.preconditioner().init(m_model->numActiveParticles(), diagonalMatrixElement, (void*)this);

	m_solver.setTolerance(m_maxError);
	m_solver.setMaxIterations(m_maxIter);
	m_solver.compute(A);

	VectorXr b(3*numParticles);
	VectorXr x(3*numParticles);
	VectorXr g(3*numParticles);

    computeRHS(b, g);

    //////////////////////////////////////////////////////////////////////////
	// Solve linear system 
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("CG solve");
	x = m_solver.solveWithGuess(b, g);
	m_iterations = (int)m_solver.iterations();
	STOP_TIMING_AVG;
	INCREASE_COUNTER("Visco iterations", static_cast<Real>(m_iterations));

    applyForces(x);
}

void Viscosity_Casson::applyForces(const VectorXr &x) 
{
    const int numParticles = (int) m_model->numActiveParticles();
    Simulation* sim = Simulation::getCurrent();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    const Real h = sim->getSupportRadius();
    const Real h2 = h*h;
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();
    const Real density0 = m_model->getDensity0();
    const Real mu = m_viscosity * density0;
    const Real mub = m_boundaryViscosity * density0;
    const Real sphereVolume = static_cast<Real>(4.0 / 3.0 * M_PI) * h2*h;
    Real d = 10.0;
    if (sim->is2DSimulation())
        d = 8.0;

	#pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static)
        for (int i = 0; i < (int)numParticles; i++)
        {
            // Compute the acceleration from the velocity change
            Vector3r &ai = m_model->getAcceleration(i);
            const Vector3r newVi(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
            ai += (1.0 / dt) * (newVi - m_model->getVelocity(i));
            m_vDiff[i] = (newVi - m_model->getVelocity(i));

            const Vector3r &xi = m_model->getPosition(i);
            const Real density_i = m_model->getDensity(i);
            const Real m_i = m_model->getMass(i);

            //////////////////////////////////////////////////////////////////////////
            // Boundary
            //////////////////////////////////////////////////////////////////////////
            if (mub != 0.0)
            {
                if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
                {
                    forall_boundary_neighbors(
                            const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
                            const Vector3r xixj = xi - xj;
                            const Vector3r gradW = sim->gradW(xixj);
                            const Vector3r a = d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (newVi - vj).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
                            bm_neighbor->addForce(xj, -m_model->getMass(i) / density_i * a);
                    );
                }
                else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
                {
                    forall_density_maps(
                            const Vector3r xixj = xi - xj;
                            Vector3r normal = -xixj;
                            const Real nl = normal.norm();
                            if (nl > static_cast<Real>(0.0001))
                            {
                                normal /= nl;

                                Vector3r t1;
                                Vector3r t2;
                                MathFunctions::getOrthogonalVectors(normal, t1, t2);

                                const Real dist = m_tangentialDistanceFactor * h;
                                const Vector3r x1 = xj - t1 * dist;
                                const Vector3r x2 = xj + t1 * dist;
                                const Vector3r x3 = xj - t2 * dist;
                                const Vector3r x4 = xj + t2 * dist;

                                const Vector3r xix1 = xi - x1;
                                const Vector3r xix2 = xi - x2;
                                const Vector3r xix3 = xi - x3;
                                const Vector3r xix4 = xi - x4;

                                const Vector3r gradW1 = sim->gradW(xix1);
                                const Vector3r gradW2 = sim->gradW(xix2);
                                const Vector3r gradW3 = sim->gradW(xix3);
                                const Vector3r gradW4 = sim->gradW(xix4);

                                // each sample point represents the quarter of the volume inside of the boundary
                                const Real vol = static_cast<Real>(0.25) * rho * sphereVolume;

                                Vector3r v1;
                                Vector3r v2;
                                Vector3r v3;
                                Vector3r v4;
                                bm_neighbor->getPointVelocity(x1, v1);
                                bm_neighbor->getPointVelocity(x2, v2);
                                bm_neighbor->getPointVelocity(x3, v3);
                                bm_neighbor->getPointVelocity(x4, v4);

                                // compute forces for both sample point
                                const Vector3r a1 = d * mub * vol * (newVi - v1).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
                                const Vector3r a2 = d * mub * vol * (newVi - v2).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
                                const Vector3r a3 = d * mub * vol * (newVi - v3).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
                                const Vector3r a4 = d * mub * vol * (newVi - v4).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;

                                bm_neighbor->addForce(x1, -m_model->getMass(i)/density_i * a1);
                                bm_neighbor->addForce(x2, -m_model->getMass(i)/density_i * a2);
                                bm_neighbor->addForce(x3, -m_model->getMass(i)/density_i * a3);
                                bm_neighbor->addForce(x4, -m_model->getMass(i)/density_i * a4);
                            }
                    );
                }
                else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
                {
                    forall_volume_maps(
                            const Vector3r xixj = xi - xj;
                            Vector3r normal = -xixj;
                            const Real nl = normal.norm();
                            if (nl > static_cast<Real>(0.0001))
                            {
                                normal /= nl;

                                Vector3r t1;
                                Vector3r t2;
                                MathFunctions::getOrthogonalVectors(normal, t1, t2);

                                const Real dist = m_tangentialDistanceFactor * h;
                                const Vector3r x1 = xj - t1*dist;
                                const Vector3r x2 = xj + t1*dist;
                                const Vector3r x3 = xj - t2*dist;
                                const Vector3r x4 = xj + t2*dist;

                                const Vector3r xix1 = xi - x1;
                                const Vector3r xix2 = xi - x2;
                                const Vector3r xix3 = xi - x3;
                                const Vector3r xix4 = xi - x4;

                                const Vector3r gradW1 = sim->gradW(xix1);
                                const Vector3r gradW2 = sim->gradW(xix2);
                                const Vector3r gradW3 = sim->gradW(xix3);
                                const Vector3r gradW4 = sim->gradW(xix4);

                                // each sample point represents the quarter of the volume inside of the boundary
                                const Real vol = static_cast<Real>(0.25) * Vj;

                                Vector3r v1;
                                Vector3r v2;
                                Vector3r v3;
                                Vector3r v4;
                                bm_neighbor->getPointVelocity(x1, v1);
                                bm_neighbor->getPointVelocity(x2, v2);
                                bm_neighbor->getPointVelocity(x3, v3);
                                bm_neighbor->getPointVelocity(x4, v4);

                                // compute forces for both sample point
                                const Vector3r a1 = d * mub * vol * (newVi - v1).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
                                const Vector3r a2 = d * mub * vol * (newVi - v2).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
                                const Vector3r a3 = d * mub * vol * (newVi - v3).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
                                const Vector3r a4 = d * mub * vol * (newVi - v4).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;

                                bm_neighbor->addForce(x1, -m_model->getMass(i)/density_i * a1);
                                bm_neighbor->addForce(x2, -m_model->getMass(i)/density_i * a2);
                                bm_neighbor->addForce(x3, -m_model->getMass(i)/density_i * a3);
                                bm_neighbor->addForce(x4, -m_model->getMass(i)/density_i * a4);
                            }
                    );
                }
            }
        }
    }
}

void Viscosity_Casson::computeRHS(VectorXr &b, VectorXr &g) 
{

    const int numParticles = (int) m_model->numActiveParticles();
    Simulation* sim = Simulation::getCurrent();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    const Real h = sim->getSupportRadius();
    const Real h2 = h*h;
    const Real dt = TimeManager::getCurrent()->getTimeStepSize();
    const Real density0 = m_model->getDensity0();
    const Real mu = m_viscosity * density0;
    const Real mub = m_boundaryViscosity * density0;
    const Real sphereVolume = static_cast<Real>(4.0 / 3.0 * M_PI) * h2*h;
    Real d = 10.0;
    if (sim->is2DSimulation())
        d = 8.0;

    //////////////////////////////////////////////////////////////////////////
    // Compute RHS
    //////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
    {
        #pragma omp for schedule(static) nowait
        for (int i = 0; i < (int)numParticles; i++)
        {
            const Vector3r &vi = m_model->getVelocity(i);
            const Vector3r &xi = m_model->getPosition(i);
            const Real density_i = m_model->getDensity(i);
            const Real m_i = m_model->getMass(i);
            Vector3r bi = Vector3r::Zero();


            //////////////////////////////////////////////////////////////////////////
            // Boundary
            //////////////////////////////////////////////////////////////////////////
            if (mub != 0.0)
            {
                if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
                {
                    forall_boundary_neighbors(
                            const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
                            const Vector3r xixj = xi - xj;
                            const Vector3r gradW = sim->gradW(xixj);
                            const Vector3r a = d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vj).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
                            bi += a;
                    );
                }
                else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
                {
                    forall_density_maps(
                            const Vector3r xixj = xi - xj;
                            Vector3r normal = -xixj;
                            const Real nl = normal.norm();
                            if (nl > static_cast<Real>(0.0001))
                            {
                                normal /= nl;

                                Vector3r t1;
                                Vector3r t2;
                                MathFunctions::getOrthogonalVectors(normal, t1, t2);

                                const Real dist = m_tangentialDistanceFactor * h;
                                const Vector3r x1 = xj - t1 * dist;
                                const Vector3r x2 = xj + t1 * dist;
                                const Vector3r x3 = xj - t2 * dist;
                                const Vector3r x4 = xj + t2 * dist;

                                const Vector3r xix1 = xi - x1;
                                const Vector3r xix2 = xi - x2;
                                const Vector3r xix3 = xi - x3;
                                const Vector3r xix4 = xi - x4;

                                const Vector3r gradW1 = sim->gradW(xix1);
                                const Vector3r gradW2 = sim->gradW(xix2);
                                const Vector3r gradW3 = sim->gradW(xix3);
                                const Vector3r gradW4 = sim->gradW(xix4);

                                // each sample point represents the quarter of the volume inside of the boundary
                                const Real vol = static_cast<Real>(0.25) * rho * sphereVolume;

                                Vector3r v1;
                                Vector3r v2;
                                Vector3r v3;
                                Vector3r v4;
                                bm_neighbor->getPointVelocity(x1, v1);
                                bm_neighbor->getPointVelocity(x2, v2);
                                bm_neighbor->getPointVelocity(x3, v3);
                                bm_neighbor->getPointVelocity(x4, v4);

                                // compute forces for both sample point
                                const Vector3r a1 = d * mub * vol * (v1).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
                                const Vector3r a2 = d * mub * vol * (v2).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
                                const Vector3r a3 = d * mub * vol * (v3).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
                                const Vector3r a4 = d * mub * vol * (v4).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
                                bi += a1 + a2 + a3 + a4;
                            }
                    );
                }
                else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
                {
                    forall_volume_maps(
                            const Vector3r xixj = xi - xj;
                            Vector3r normal = -xixj;
                            const Real nl = normal.norm();
                            if (nl > static_cast<Real>(0.0001))
                            {
                                normal /= nl;

                                Vector3r t1;
                                Vector3r t2;
                                MathFunctions::getOrthogonalVectors(normal, t1, t2);

                                const Real dist = m_tangentialDistanceFactor * h;
                                const Vector3r x1 = xj - t1*dist;
                                const Vector3r x2 = xj + t1*dist;
                                const Vector3r x3 = xj - t2*dist;
                                const Vector3r x4 = xj + t2*dist;

                                const Vector3r xix1 = xi - x1;
                                const Vector3r xix2 = xi - x2;
                                const Vector3r xix3 = xi - x3;
                                const Vector3r xix4 = xi - x4;

                                const Vector3r gradW1 = sim->gradW(xix1);
                                const Vector3r gradW2 = sim->gradW(xix2);
                                const Vector3r gradW3 = sim->gradW(xix3);
                                const Vector3r gradW4 = sim->gradW(xix4);

                                // each sample point represents the quarter of the volume inside of the boundary
                                const Real vol = static_cast<Real>(0.25) * Vj;

                                Vector3r v1;
                                Vector3r v2;
                                Vector3r v3;
                                Vector3r v4;
                                bm_neighbor->getPointVelocity(x1, v1);
                                bm_neighbor->getPointVelocity(x2, v2);
                                bm_neighbor->getPointVelocity(x3, v3);
                                bm_neighbor->getPointVelocity(x4, v4);

                                // compute forces for both sample point
                                const Vector3r a1 = d * mub * vol * (v1).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
                                const Vector3r a2 = d * mub * vol * (v2).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
                                const Vector3r a3 = d * mub * vol * (v3).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
                                const Vector3r a4 = d * mub * vol * (v4).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
                                bi += a1 + a2 + a3 + a4;
                            }
                    );
                }
            }

            b[3*i] = vi[0] - dt/density_i * bi[0];
            b[3*i+1] = vi[1] - dt/density_i * bi[1];
            b[3*i+2] = vi[2] - dt/density_i * bi[2];

            // Warmstart
            g[3 * i] = vi[0] + m_vDiff[i][0];
            g[3 * i + 1] = vi[1] + m_vDiff[i][1];
            g[3 * i + 2] = vi[2] + m_vDiff[i][2];
        }
    }
}


void Viscosity_Casson::reset()
{
}

void Viscosity_Casson::performNeighborhoodSearchSort()
{
}

