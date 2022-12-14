target_sources(SPlisHSPlasH 
PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/SurfaceParticles/SurfaceParticles.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SurfaceParticles/SurfaceParticles.h

    ${CMAKE_CURRENT_LIST_DIR}/Viscosity_Casson/Viscosity_Casson.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Viscosity_Casson/Viscosity_Casson.h

    ${CMAKE_CURRENT_LIST_DIR}/Diffusion/Coagulation.h	
    ${CMAKE_CURRENT_LIST_DIR}/Diffusion/Coagulation.cpp

	${CMAKE_CURRENT_LIST_DIR}/Diffusion/TemperatureDiffusion.h
    ${CMAKE_CURRENT_LIST_DIR}/Diffusion/TemperatureDiffusion.cpp

    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton.h
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton.cpp

    ${CMAKE_CURRENT_LIST_DIR}/MyPartioReader/MyPartioReader.cpp 
	${CMAKE_CURRENT_LIST_DIR}/MyPartioReader/MyPartioReader.h

    ${CMAKE_CURRENT_LIST_DIR}/ShapeMatching/ShapeMatching.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ShapeMatching/ShapeMatching.h
)




# set(PLASTIC_HEADER_FILES	
# 	My/Plastic/Plastic.h	#add plastic
# 	)

# set(PLASTIC_SOURCE_FILES	
# 	My/Plastic/Plastic.cpp	#add plastic
# 	)