target_sources(SPlisHSPlasH 
PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/findSurfaceParticles/SurfaceParticles.cpp

    ${CMAKE_CURRENT_LIST_DIR}/Viscosity_Casson/Viscosity_Casson.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Viscosity_Casson/Viscosity_Casson.h
    ${CMAKE_CURRENT_LIST_DIR}/Diffusion/Coagulation.h	
    ${CMAKE_CURRENT_LIST_DIR}/Diffusion/Coagulation.cpp
	${CMAKE_CURRENT_LIST_DIR}/Diffusion/TemperatureDiffusion.h
    ${CMAKE_CURRENT_LIST_DIR}/Diffusion/TemperatureDiffusion.cpp
)


set(DIFFUSION_HEADER_FILES	#add diffusion
	Diffusion/Coagulation.h	
	Diffusion/TemperatureDiffusion.h
	)
	
set(DIFFUSION_SOURCE_FILES	#add diffusion
	Diffusion/Coagulation.cpp
	Diffusion/TemperatureDiffusion.cpp
	)


# set(NON_NEWTON_HEADER_FILES
# 	NonNewton/NonNewton.h	#add NonNewton
# 	)
	
# set(NON_NEWTON_SOURCE_FILES
# 	NonNewton/NonNewton.cpp	#add NonNewton
# 	)


# set(RIGID_BODY_SOURCE_FILES	#add rigid body
# 	RigidBody/RigidBody.cpp
# 	)
	
# set(RIGID_BODY_HEADER_FILES	
# 	RigidBody/RigidBody.h	#add rigid body
# 	)

# set(PLASTIC_HEADER_FILES	
# 	My/Plastic/Plastic.h	#add plastic
# 	)

# set(PLASTIC_SOURCE_FILES	
# 	My/Plastic/Plastic.cpp	#add plastic
# 	)