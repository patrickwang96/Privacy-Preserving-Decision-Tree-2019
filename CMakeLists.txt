cmake_minimum_required(VERSION 3.10)
project(Privacy_Preserving_Decision_Tree_2019)

set(CMAKE_CXX_STANDARD 11)

find_library(gmp required)
include_directories(/usr/include/eigen3/ /usr/include/ /usr/local/include/ ../libOTe/ ../libOTe/cryptoTools ../libOTe/cryptoTools/thirdparty/linux/boost ../libOTe/cryptoTools/thirdparty/linux/miracl)
link_directories(/usr/local/lib/ ../libOTe/lib/ ../libOTe/cryptoTools/thirdparty/linux/boost/stage/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "./build")

set(COMMON_SRCS types.cc types.h config.cc config.h secret_sharing.cc secret_sharing.h utils.cc utils.h)
set(CLIENT_SRCS test_cloud_client.cpp test_cloud_client.h client_main.cpp client_decision_tree.cpp client_decision_tree.h)
set(SERVER_SRCS test_cloud_server.cpp test_cloud_server.h server_main.cpp server_decision_tree.cpp server_decision_tree.h)
set(MAIN_SRCS test_driver.h test_driver.cc decision_tree.h decision_tree.cc main.cc)
set(NET_SRCS network.cc network.h secret_sharing_efficient_tools.h secret_sharing_efficient_tools.cpp bit_compression.cpp)


add_executable(main ${COMMON_SRCS} ${MAIN_SRCS})
add_executable(client ${COMMON_SRCS} ${CLIENT_SRCS} ${NET_SRCS})
add_executable(server ${COMMON_SRCS} ${SERVER_SRCS} ${NET_SRCS})
target_link_libraries(main gmp gmpxx)
target_link_libraries(client gmp gmpxx)
target_link_libraries(server gmp gmpxx)
