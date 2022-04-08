
if(NOT "D:/codes/SPH/SPlisHSPlasH/SPlisHSPlasH/extern/cuNSearch/src/Ext_NeighborhoodSearch-stamp/Ext_NeighborhoodSearch-gitinfo.txt" IS_NEWER_THAN "D:/codes/SPH/SPlisHSPlasH/SPlisHSPlasH/extern/cuNSearch/src/Ext_NeighborhoodSearch-stamp/Ext_NeighborhoodSearch-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: 'D:/codes/SPH/SPlisHSPlasH/SPlisHSPlasH/extern/cuNSearch/src/Ext_NeighborhoodSearch-stamp/Ext_NeighborhoodSearch-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "D:/codes/SPH/SPlisHSPlasH/SPlisHSPlasH/extern/cuNSearch/src/Ext_NeighborhoodSearch"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: 'D:/codes/SPH/SPlisHSPlasH/SPlisHSPlasH/extern/cuNSearch/src/Ext_NeighborhoodSearch'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "D:/App/Git/cmd/git.exe"  clone --no-checkout --config "advice.detachedHead=false" "https://github.com/InteractiveComputerGraphics/cuNSearch.git" "Ext_NeighborhoodSearch"
    WORKING_DIRECTORY "D:/codes/SPH/SPlisHSPlasH/SPlisHSPlasH/extern/cuNSearch/src"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/InteractiveComputerGraphics/cuNSearch.git'")
endif()

execute_process(
  COMMAND "D:/App/Git/cmd/git.exe"  checkout aba3da18cb4f45cd05d729465d1725891ffc33da --
  WORKING_DIRECTORY "D:/codes/SPH/SPlisHSPlasH/SPlisHSPlasH/extern/cuNSearch/src/Ext_NeighborhoodSearch"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'aba3da18cb4f45cd05d729465d1725891ffc33da'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "D:/App/Git/cmd/git.exe"  submodule update --recursive --init 
    WORKING_DIRECTORY "D:/codes/SPH/SPlisHSPlasH/SPlisHSPlasH/extern/cuNSearch/src/Ext_NeighborhoodSearch"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: 'D:/codes/SPH/SPlisHSPlasH/SPlisHSPlasH/extern/cuNSearch/src/Ext_NeighborhoodSearch'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "D:/codes/SPH/SPlisHSPlasH/SPlisHSPlasH/extern/cuNSearch/src/Ext_NeighborhoodSearch-stamp/Ext_NeighborhoodSearch-gitinfo.txt"
    "D:/codes/SPH/SPlisHSPlasH/SPlisHSPlasH/extern/cuNSearch/src/Ext_NeighborhoodSearch-stamp/Ext_NeighborhoodSearch-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: 'D:/codes/SPH/SPlisHSPlasH/SPlisHSPlasH/extern/cuNSearch/src/Ext_NeighborhoodSearch-stamp/Ext_NeighborhoodSearch-gitclone-lastrun.txt'")
endif()

