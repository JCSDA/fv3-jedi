list( APPEND fv3jedi_test_data
Data/C96_RESTART_2016-01-01-06z/C96_mosaic.nc
Data/C96_RESTART_2016-01-01-06z/coupler.res
Data/C96_RESTART_2016-01-01-06z/fv_core.res.nc
Data/C96_RESTART_2016-01-01-06z/grid_spec.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/C96_grid_spec.tile1.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/C96_grid_spec.tile2.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/C96_grid_spec.tile3.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/C96_grid_spec.tile4.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/C96_grid_spec.tile5.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/C96_grid_spec.tile6.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/C96_grid.tile1.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/C96_grid.tile2.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/C96_grid.tile3.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/C96_grid.tile4.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/C96_grid.tile5.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/C96_grid.tile6.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/fv_core.res.tile1.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/fv_core.res.tile2.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/fv_core.res.tile3.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/fv_core.res.tile4.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/fv_core.res.tile5.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/fv_core.res.tile6.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/fv_tracer.res.tile1.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/fv_tracer.res.tile2.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/fv_tracer.res.tile3.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/fv_tracer.res.tile4.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/fv_tracer.res.tile5.nc
Data/C96_RESTART_2016-01-01-06z/INPUT/fv_tracer.res.tile6.nc
)

list( APPEND fv3jedi_test_fv3files
Data/fv3files/field_table
Data/fv3files/fmsmpp.nml
Data/fv3files/input.nml
Data/fv3files/input_c48.nml
Data/fv3files/inputpert_fsoi.nml
Data/fv3files/inputpert_4dvar.nml
)

list( APPEND fv3jedi_test_restart
Data/RESTART/PlotCubeSphereField.py
)


list( APPEND fv3jedi_test_input
testinput/geometry.json
testinput/state.json
testinput/analytic_init.json
testinput/model.json
testinput/linearmodel.json
testinput/3dvar.json
testinput/4dvar.drplanczos.json
testinput/forecast.json
testinput/hofx.json
testinput/dirac.nicas.json
testinput/analytic_init.json
testinput/increment.json
testinput/4denvar.nicas.json
testinput/errorcov.json
testinput/localization.json
compare.sh
)

list( APPEND fv3jedi_testoutput
testoutput/hofx.test  
testoutput/3dvar.test
testoutput/4dvar.drplanczos.test
testoutput/dirac.nicas.test
)

# Create Data directory for model input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/C96_RESTART_2016-01-01-06z)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/C96_RESTART_2016-01-01-06z/INPUT)
foreach(FILENAME ${fv3jedi_test_data})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for model input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/fv3files)
foreach(FILENAME ${fv3jedi_test_fv3files})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for model input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/RESTART)
foreach(FILENAME ${fv3jedi_test_restart})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${fv3jedi_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for reference output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
foreach(FILENAME ${fv3jedi_testoutput})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

ecbuild_add_resources( TARGET   fv3jedi_test_scripts
                       SOURCES_PACK
                       ${fv3jedi_test_input}
                     )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

#####################################################################

# Test interface classes with FV3JEDI

ecbuild_add_test( TARGET  test_fv3jedi_geometry
                  BOOST
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/geometry.json"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_state
                  BOOST
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/state.json"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_aninit
                  BOOST
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/analytic_init.json"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_model
                  BOOST
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model.json"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_linearmodel
                  BOOST
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/linearmodel.json"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_modelaux
                  BOOST
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    "testinput/model.json"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_increment
                  BOOST
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment.json"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_errorcovariance
                  BOOST
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/errorcov.json"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_modelauxcovariance
                  BOOST
                  SOURCES executables/TestModelAuxCovariance.cc
                  ARGS    "testinput/model.json"
                  MPI     6
                  LIBS    fv3jedi )

ecbuild_add_test( TARGET  test_fv3jedi_localization
                  BOOST
                  SOURCES executables/TestLocalization.cc
                  ARGS    "testinput/localization.json"
                  MPI     6
                  LIBS    fv3jedi )

# Test oops applications with FV3JEDI
# compare.sh does not currently work with mpirun

#ecbuild_add_test( TARGET test_fv3jedi_dirac_nicas
#                  TYPE SCRIPT
#                  MPI     6
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/fv3jedi_dirac.x testinput/dirac.nicas.json"
#                       testoutput/dirac.nicas.test
#                  DEPENDS fv3jedi_dirac.x )
#
#ecbuild_add_test( TARGET test_fv3jedi_hofx
#                  TYPE SCRIPT
#                  MPI     6
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/fv3jedi_hofx.x testinput/hofx.json"
#                       testoutput/hofx.test
#                  DEPENDS fv3jedi_hofx.x )
#
#ecbuild_add_test( TARGET test_fv3jedi_3dvar
#                  TYPE SCRIPT
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/fv3jedi_4dvar.x testinput/3dvar.json"
#                       testoutput/3dvar.test
#                  DEPENDS fv3jedi_4dvar.x )
#
#ecbuild_add_test( TARGET test_fv3jedi_3denvar
#                  TYPE SCRIPT
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/fv3jedi_4dvar.x testinput/3dvar.json"
#                       testoutput/3denvar_nicas.test
#                  DEPENDS fv3jedi_4dvar.x )
#
#ecbuild_add_test( TARGET test_fv3jedi_4dvar_drplanczos
#                  TYPE SCRIPT
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/fv3jedi_4dvar.x testinput/4dvar.drplanczos.json"
#                       testoutput/4dvar.drplanczos.test
#                  DEPENDS fv3jedi_4dvar.x )
