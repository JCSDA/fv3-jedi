# (C) Copyright 2018-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Extra macros to eliminate repetition

# Macro to link list of files from source to destination
macro( LINK_FILES filelist src_dir dst_dir )
  foreach(FILENAME ${filelist})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${src_dir}/${FILENAME}
      ${dst_dir}/${FILENAME}
    )
  endforeach(FILENAME)
endmacro()

macro( LINK_FILES_DIR filelist dst_dir )
  foreach(FILENAME ${filelist})
    execute_process( COMMAND ln -sf ${FILENAME} ${dst_dir} )
  endforeach(FILENAME)
endmacro()

# macro to create a symlink from src to dst
function(CREATE_SYMLINK src dst)
    foreach (FILENAME ${ARGN})
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${src}/${FILENAME}
            ${dst}/${FILENAME} )
        endforeach(FILENAME)
endfunction(CREATE_SYMLINK)

# macro to create a symlink from src to dst with just filename
function(CREATE_SYMLINK_FILENAME src dst)
    foreach (FILENAME ${ARGN})
        get_filename_component(filename ${FILENAME} NAME )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${src}/${FILENAME}
            ${dst}/${filename} )
        endforeach(FILENAME)
endfunction(CREATE_SYMLINK_FILENAME)

# Macro to include GFS restart files in testing
macro( GFS_FILES_BKG path date )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.coupler.res )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile1.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile2.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile3.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile4.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile5.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile6.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_srf_wnd.res.tile1.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_srf_wnd.res.tile2.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_srf_wnd.res.tile3.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_srf_wnd.res.tile4.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_srf_wnd.res.tile5.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_srf_wnd.res.tile6.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile1.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile2.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile3.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile4.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile5.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile6.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.sfc_data.tile1.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.sfc_data.tile2.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.sfc_data.tile3.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.sfc_data.tile4.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.sfc_data.tile5.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.sfc_data.tile6.nc )
endmacro()

macro( GFS_FILES_ENS path date )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.coupler.res )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile1.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile2.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile3.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile4.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile5.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile6.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile1.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile2.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile3.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile4.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile5.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile6.nc )
endmacro()

macro( GFS_AERO_FILES_BKG path date )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.coupler.res )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile1.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile2.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile3.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile4.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile5.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_core.res.tile6.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile1.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile2.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile3.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile4.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile5.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile6.nc )
endmacro()

macro( GFS_AERO_FILES_ENS path date )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.coupler.res )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile1.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile2.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile3.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile4.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile5.nc )
    list( APPEND fv3jedi_gfs_test_data ${path}/${date}.fv_tracer.res.tile6.nc )
endmacro()

macro( LAM_CMAQ_FILES_BKG path date)
    list( APPEND fv3jedi_lam_test_data ${path}/${date}.coupler.res )
    list( APPEND fv3jedi_lam_test_data ${path}/${date}.fv_core.res.tile1.nc )
    list( APPEND fv3jedi_lam_test_data ${path}/${date}.fv_srf_wnd.res.tile1.nc )
    list( APPEND fv3jedi_lam_test_data ${path}/${date}.fv_tracer.res.tile1.nc )
    list( APPEND fv3jedi_lam_test_data ${path}/${date}.sfc_data.nc )
endmacro()
