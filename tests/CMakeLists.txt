add_executable(testGlobalElement
    mesh/testGlobalElement.cpp
)

target_link_libraries(testGlobalElement
    PRIVATE globalElement
)
