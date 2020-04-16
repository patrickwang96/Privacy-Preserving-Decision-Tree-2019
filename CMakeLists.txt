cmake_minimum_required(VERSION 3.16)
project(Privacy_Preserving_Decision_Tree_2019)

set(CMAKE_CXX_STANDARD 11)

find_library(gmp required)
include_directories(/usr/local/include/eigen3/ /usr/local/include/ ../libOTe/ ../libOTe/cryptoTools ../libOTe/cryptoTools/thirdparty/linux/boost ../libOTe/cryptoTools/thirdparty/linux/miracl)
link_directories(/usr/local/lib/ ../libOTe/lib/ ../libOTe/cryptoTools/thirdparty/linux/boost/stage/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "./build")

set(COMMON_SRCS types.cc types.h config.cc config.h secret_sharing.cc secret_sharing.h utils.cc utils.h test_driver.cc test_driver.h decision_tree.cc decision_tree.h)


add_executable(main ${COMMON_SRCS} main.cc)
target_link_libraries(main gmp gmpxx)
