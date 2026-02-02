ecbuild_add_executable( TARGET  fv3jedi_addincrement.x
                        SOURCES fv3jediAddIncrement.cc
                        LIBS	${FV3JEDI_LIBRARIES}
                      )

ecbuild_add_executable( TARGET  fv3jedi_convertincrement.x
                        SOURCES fv3jediConvertIncrement.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                                saber
                       )

ecbuild_add_executable( TARGET  fv3jedi_controlpert.x
                        SOURCES fv3jediControlPert.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                                saber
                       )

ecbuild_add_executable( TARGET  fv3jedi_convertstate.x
                        SOURCES fv3jediConvertState.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                                saber
                       )

ecbuild_add_executable( TARGET  fv3jedi_converttostructuredgrid.x
                        SOURCES fv3jediConvertToStructuredGrid.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                       )

ecbuild_add_executable( TARGET  fv3jedi_diffstates.x
                        SOURCES fv3jediDiffStates.cc
                        LIBS	${FV3JEDI_LIBRARIES}
                      )

ecbuild_add_executable( TARGET  fv3jedi_eda.x
                        SOURCES fv3jediEDA.cc
                        LIBS	${FV3JEDI_LIBRARIES}
                                saber
                      )

ecbuild_add_executable( TARGET  fv3jedi_enshofx.x
                        SOURCES fv3jediEnsHofX.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                      )

ecbuild_add_executable( TARGET  fv3jedi_ensrecenter.x
                        SOURCES fv3jediEnsRecenter.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                        )

ecbuild_add_executable( TARGET  fv3jedi_ensmeanvariance.x
                        SOURCES fv3jediEnsMeanVariance.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                        )

ecbuild_add_executable( TARGET  fv3jedi_forecast.x
                        SOURCES fv3jediForecast.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                       )

ecbuild_add_executable( TARGET  fv3jedi_EnsGETKF.x
                        SOURCES fv3jediEnsGETKF.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                       )

ecbuild_add_executable( TARGET  fv3jedi_gen_hybrid_linear_model_coeffs.x
                        SOURCES fv3jediGenHybridLinearModelCoeffs.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                       )

ecbuild_add_executable( TARGET  fv3jedi_hofx.x
                        SOURCES fv3jediHofX.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                      )

ecbuild_add_executable( TARGET  fv3jedi_hofx_nomodel.x
                        SOURCES fv3jediHofXNoModel.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                      )

ecbuild_add_executable( TARGET  fv3jedi_letkf.x
                        SOURCES fv3jediLETKF.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                      )

ecbuild_add_executable( TARGET  fv3jedi_error_covariance_toolbox.x
                        SOURCES fv3jediErrorCovarianceToolbox.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                                saber
                      )

ecbuild_add_executable( TARGET  fv3jedi_var.x
                        SOURCES fv3jediVar.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                                saber
                      )

ecbuild_add_executable( TARGET  fv3jedi_adjointforecast.x
                        SOURCES fv3jediAdjointForecast.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                      )

ecbuild_add_executable( TARGET  fv3jedi_genenspertB.x
                        SOURCES fv3jediGenEnsPertB.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                                saber
                      )

ecbuild_add_executable( TARGET  fv3jedi_processperts.x
                        SOURCES fv3jediProcessPerts.cc
                        LIBS    ${FV3JEDI_LIBRARIES}
                                saber
                      )

ecbuild_add_executable( TARGET  fv3jedi_cold_to_warm.x
                        SOURCES fv3jediColdToWarm.f90
                        LIBS	${FV3JEDI_LIBRARIES}
                      )

ecbuild_add_executable( TARGET  fv3jedi_a_to_d.x
                        SOURCES fv3jediAtoD.f90
                        LIBS	${FV3JEDI_LIBRARIES}
                      )

ecbuild_add_executable( TARGET  fv3jedi_d_to_a.x
                        SOURCES fv3jediDtoA.f90
                        LIBS	${FV3JEDI_LIBRARIES}
                      )

ecbuild_add_executable( TARGET  fv3jedi_geos_to_fms.x
                        SOURCES fv3jediGeosToFms.f90
                        LIBS	${FV3JEDI_LIBRARIES}
                      )
ecbuild_add_executable( TARGET fv3jedi_tlm_toolbox.x
                        SOURCES fv3jediTLMToolbox.cc
                        LIBS   ${FV3JEDI_LIBRARIES}
                      )

# Plotting utility
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fv3jedi_plot_field.py ${CMAKE_BINARY_DIR}/bin/fv3jedi_plot_field.x COPYONLY)
execute_process( COMMAND chmod +x ${CMAKE_BINARY_DIR}/bin/fv3jedi_plot_field.x)
