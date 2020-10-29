ecbuild_add_executable( TARGET  fv3jedi_addincrement.x
                        SOURCES fv3jediAddIncrement.cc
                        LIBS	fv3jedi
                      )

ecbuild_add_executable( TARGET  fv3jedi_convertstate.x
                        SOURCES fv3jediConvertState.cc
                        LIBS    fv3jedi
                                saber
                       )

ecbuild_add_executable( TARGET  fv3jedi_diffstates.x
                        SOURCES fv3jediDiffStates.cc
                        LIBS	fv3jedi
                      )

ecbuild_add_executable( TARGET  fv3jedi_dirac.x
                        SOURCES fv3jediDirac.cc
                        LIBS	fv3jedi
                                saber
                      )

ecbuild_add_executable( TARGET  fv3jedi_eda.x
                        SOURCES fv3jediEDA.cc
                        LIBS	fv3jedi
                                saber
                      )

ecbuild_add_executable( TARGET  fv3jedi_enshofx.x
                        SOURCES fv3jediEnsHofX.cc
                        LIBS    fv3jedi
                      )

ecbuild_add_executable( TARGET  fv3jedi_ensvariance.x
                        SOURCES fv3jediEnsVariance.cc
                        LIBS    fv3jedi
                        )

ecbuild_add_executable( TARGET  fv3jedi_forecast.x
                        SOURCES fv3jediForecast.cc
                        LIBS    fv3jedi
                       )

ecbuild_add_executable( TARGET  fv3jedi_hofx.x
                        SOURCES fv3jediHofX.cc
                        LIBS    fv3jedi
                      )

ecbuild_add_executable( TARGET  fv3jedi_hofx_nomodel.x
                        SOURCES fv3jediHofXNoModel.cc
                        LIBS    fv3jedi
                      )

ecbuild_add_executable( TARGET  fv3jedi_letkf.x
                        SOURCES fv3jediLETKF.cc
                        LIBS    fv3jedi
                      )

ecbuild_add_executable( TARGET  fv3jedi_parameters.x
                        SOURCES fv3jediEstimateParams.cc
                        LIBS    fv3jedi
                                saber
                      )

ecbuild_add_executable( TARGET  fv3jedi_var.x
                        SOURCES fv3jediVar.cc
                        LIBS    fv3jedi
                                saber
                      )
#ecbuild_add_executable( TARGET  fv3jedi_linearforecast.x
#                        SOURCES fv3jediLinearForecast.cc
#                        LIBS    fv3jedi
#                       )

# Plotting utility
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fv3jedi_plot_field.py ${CMAKE_BINARY_DIR}/bin/fv3jedi_plot_field.x COPYONLY)
execute_process( COMMAND chmod +x ${CMAKE_BINARY_DIR}/bin/fv3jedi_plot_field.x)
