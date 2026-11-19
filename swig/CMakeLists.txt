include(OpenEyeSWIG)

# Ask SWIG for the headers the .i file includes, so editing a public header
# regenerates the wrapper (CMake 3.20+; ignored by older versions).
set(SWIG_USE_SWIG_DEPENDENCIES TRUE)

openeye_add_swig_module(
    NAME oehealth
    SWIG_FILE oehealth.i
    LINK_LIBS oehealth
    PYTHON_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/python/oehealth
    STABLE_ABI ${OEHEALTH_USE_STABLE_ABI}
    SWIG_FLAGS "-doxygen" "-I${CMAKE_SOURCE_DIR}/include"
    EXPECTED_LIB_VARS
        OECHEM_LIBRARY OESYSTEM_LIBRARY OEPLATFORM_LIBRARY
        OEMATH_LIBRARY OEZSTD_LIBRARY OEBIO_LIBRARY OEGRID_LIBRARY
    COMPILE_DEFS $<$<PLATFORM_ID:Windows>:NOGDI>
    INIT_PY ${CMAKE_SOURCE_DIR}/python/oehealth/__init__.py
)
