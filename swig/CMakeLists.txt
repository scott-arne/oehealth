include(OpenEyeSWIG)

openeye_add_swig_module(
    NAME oehealth
    SWIG_FILE oehealth.i
    LINK_LIBS oehealth
    PYTHON_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/python/oehealth
    STABLE_ABI ${OEHEALTH_USE_STABLE_ABI}
    SWIG_FLAGS "-doxygen"
    EXPECTED_LIB_VARS
        OECHEM_LIBRARY OESYSTEM_LIBRARY OEPLATFORM_LIBRARY
        OEMATH_LIBRARY OEZSTD_LIBRARY
    INIT_PY ${CMAKE_SOURCE_DIR}/python/oehealth/__init__.py
)
