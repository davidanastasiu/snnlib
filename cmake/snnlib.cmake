include_guard()

# get library version from code
function(snnlib_getversion version_arg soversion_arg)
    # Parse the current version from the snnlib header
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/snnlib/version.h" snnlib_version_defines
            REGEX "#define SNNLIB_VERSION_(MAJOR|MINOR|PATCH)")
    foreach(ver ${snnlib_version_defines})
        if(ver MATCHES "#define SNNLIB_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
            set(SNNLIB_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
        endif()
    endforeach()
    set(VERSION ${SNNLIB_VERSION_MAJOR}.${SNNLIB_VERSION_MINOR}.${SNNLIB_VERSION_PATCH})

    message(STATUS "snnlib version ${VERSION}")

    # Return the information to the caller
    set(${version_arg} ${VERSION} PARENT_SCOPE)
    set(${soversion_arg} ${SNNLIB_VERSION_MAJOR} PARENT_SCOPE)
endfunction()