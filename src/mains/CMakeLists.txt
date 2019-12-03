ecbuild_add_executable( TARGET  fv3jedi_addincrement.x
                        SOURCES fv3jediAddIncrement.cc
                        LIBS	fv3jedi
                      )

ecbuild_add_executable( TARGET  fv3jedi_convertstate.x
                        SOURCES fv3jediConvertState.cc
                        LIBS    fv3jedi
                       )

ecbuild_add_executable( TARGET  fv3jedi_diffstates.x
                        SOURCES fv3jediDiffStates.cc
                        LIBS	fv3jedi
                      )

ecbuild_add_executable( TARGET  fv3jedi_dirac.x
                        SOURCES fv3jediDirac.cc
                        LIBS	fv3jedi
                      )

ecbuild_add_executable( TARGET  fv3jedi_eda.x
                        SOURCES fv3jediEDA.cc
                        LIBS	fv3jedi
                      )

ecbuild_add_executable( TARGET  fv3jedi_enshofx.x
                        SOURCES fv3jediEnsHofX.cc
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

ecbuild_add_executable( TARGET  fv3jedi_hofx3D.x
                        SOURCES fv3jediHofX3D.cc
                        LIBS    fv3jedi
                      )

ecbuild_add_executable( TARGET  fv3jedi_parameters.x
                        SOURCES fv3jediEstimateParams.cc
                        LIBS    fv3jedi
                      )

ecbuild_add_executable( TARGET  fv3jedi_var.x
                        SOURCES fv3jediVar.cc
                        LIBS    fv3jedi
                      )

#ecbuild_add_executable( TARGET  fv3jedi_linearforecast.x
#                        SOURCES fv3jediLinearForecast.cc
#                        LIBS    fv3jedi
#                       )
