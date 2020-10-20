Project(BNLD)
cmake_minimum_required (VERSION 2.6)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has  C++11 support.")
else()
	message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

SET(SRCs
	BNLD_main
	IncrOwen_main
)

FOREACH(FILE ${SRCs})
  add_executable( ${FILE} ${FILE} )
  target_link_libraries ( ${FILE} )
  install(TARGETS ${FILE}
			RUNTIME DESTINATION bin
			LIBRARY DESTINATION lib
			ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)





