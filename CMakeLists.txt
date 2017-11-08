cmake_minimum_required(VERSION 3.8)
project (MQTenTen)

add_definitions(-D_UNICODE -DUNICODE)

find_package(CGAL REQUIRED COMPONENTS Core)
include( ${CGAL_USE_FILE} )

set(MQSDK_DIR "" CACHE PATH "MQInit.cpp")

if(EXISTS ${MQSDK_DIR}/mqsdk/MQInit.cpp)
  set(MQSDK_INCLUDE_DIR ${MQSDK_DIR}/mqsdk/)
elseif(EXISTS ${MQSDK_DIR}/MQInit.cpp)
  set(MQSDK_INCLUDE_DIR ${MQSDK_DIR})
endif()

if(NOT DEFINED MQSDK_INCLUDE_DIR)
  message(FATAL_ERROR "MQSDK_DIR not found!")
endif()

INCLUDE_DIRECTORIES( "${MQSDK_INCLUDE_DIR}" )

set(MQTENTEN_TENUCHI 0 CACHE BOOL "Enable Ten-uchi (GPL License Warning!)")

add_library(MQTenTen SHARED
  MQTenTen.cpp
  ${MQSDK_INCLUDE_DIR}/MQInit.cpp
  ${MQSDK_INCLUDE_DIR}/MQPlugin.cpp
  ${MQSDK_INCLUDE_DIR}/MQWidget.cpp
  ${MQSDK_INCLUDE_DIR}/MQBasePlugin.cpp
  ${MQSDK_INCLUDE_DIR}/MQSetting.cpp
  ${MQSDK_INCLUDE_DIR}/MQ3DLib.cpp
)
target_link_libraries( MQTenTen ${GMP_LIBRARIES} )
target_link_libraries( MQTenTen ${MPFR_LIBRARIES} )

IF(WIN32)
  set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /OPT:REF")
ENDIF()

if(CMAKE_CL_64)
  set_target_properties(MQTenTen PROPERTIES OUTPUT_NAME "MQTenTen64")
endif()

if(MQTENTEN_TENUCHI)
  add_definitions(-DENABLE_TENUCHI)
endif()
