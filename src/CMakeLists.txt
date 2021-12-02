set(SOURCE_FILES "${PROJECT_NAME}.cpp" "${PROJECT_NAME}_internal.cpp"
        "${PROJECT_NAME}_internal_G3.cpp" "${PROJECT_NAME}_internal_G4_simp.cpp" "${PROJECT_NAME}_internal_G4.cpp"
        PolyLog.cpp Nielsen.cpp Util.cpp Elliptic.cpp)
set(PUBLIC_HEADERS "../include/${PROJECT_NAME}.h")
set(PRIVATE_HEADERS "../include/${PROJECT_NAME}_internal.h")

add_library(${PROJECT_NAME} ${SOURCE_FILES} ${PUBLIC_HEADERS} ${PRIVATE_HEADERS})

set_target_properties("${PROJECT_NAME}" PROPERTIES PUBLIC_HEADER ${PUBLIC_HEADERS})

include(GNUInstallDirs)

set(LIBRARY_DIR "${CMAKE_INSTALL_LIBDIR}")

target_include_directories("${PROJECT_NAME}" INTERFACE
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )

install(TARGETS "${PROJECT_NAME}"
        EXPORT "${PROJECT_NAME}Targets"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        )

install(EXPORT "${PROJECT_NAME}Targets"
        FILE "${PROJECT_NAME}Targets.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
        )

include(CMakePackageConfigHelpers)

configure_package_config_file("${PROJECT_NAME}Config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
        PATH_VARS LIBRARY_DIR
        )

write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        COMPATIBILITY AnyNewerVersion
)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
        )
