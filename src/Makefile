#***************************************************************************/
#**    Author: Joey Azofeifa                                               */
#**    Date:   11/8/2015                                                   */
#**    Description: Makefile for Fast Interval Searcher                    */
#***************************************************************************/

CXX             = g++-5.1.0
CXXFLAGS        = -Wno-unused-variable -Wno-non-virtual-dtor -std=c++11 -fopenmp -Wno-write-strings -Wno-literal-suffix
EXEC            = ${PWD}/FIS
ARCH  			= getconf LONG_BIT

CPP_FLAGS_32  	= -D32_BIT
CPP_FLAGS_64  	= -D64_BIT

GCCVERSION 		= $(shell ${CXX} -dumpversion)
FIS: main.o read_in_parameters.o split.o load.o find_overlaps.o EM.o
	@printf "linking           : "
	@${CXX} ${CXXFLAGS} -o ${EXEC}  ${PWD}/main.o ${PWD}/read_in_parameters.o ${PWD}/split.o ${PWD}/load.o ${PWD}/find_overlaps.o ${PWD}/EM.o  
	@printf "done\n"
	@echo "========================================="
	@printf "finished :)\n"

main.o: 
	@echo "========================================="
	@echo "GCC version: "${GCCVERSION}
	@printf "main              : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/main.cpp 
	@printf "done\n"
EM.o: 
	@printf "EM                : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/EM.cpp 
	@printf "done\n"
load.o: 
	@printf "load              : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/load.cpp 
	@printf "done\n"

find_overlaps.o: 
	@printf "find_overlaps     : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/find_overlaps.cpp 
	@printf "done\n"

read_in_parameters.o: 
	@printf "read_in_parameters: "
	@${CXX} -c ${CXXFLAGS} ${PWD}/read_in_parameters.cpp 
	@printf "done\n"

split.o: 
	@printf "split             : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/split.cpp 
	@printf "done\n"



clean:
	@rm -f *.o

