add_executable(testReferenceElement
    mesh/testReferenceElement.cpp
)

target_link_libraries(testReferenceElement
    PRIVATE referenceElement
)
