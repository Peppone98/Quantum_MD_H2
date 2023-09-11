H2: matrix_elements.o solving_Roothan.o filling_matrices.o energy_orbitals.o Cs_evol.o der_matrix_elements.o nuclear_matrices.o X_evol.o CG_routine.o DFT_xc.o DFT_derivatives.o Adaptive.o file_handling.o main_H2.o
	c++ -o H2 matrix_elements.o solving_Roothan.o filling_matrices.o energy_orbitals.o Cs_evol.o der_matrix_elements.o nuclear_matrices.o X_evol.o CG_routine.o DFT_xc.o DFT_derivatives.o Adaptive.o file_handling.o main_H2.o -lm -lgsl -lgslcblas -lxc

matrix_elements.o: matrix_elements.cpp definitions.h
	c++ -c matrix_elements.cpp

solving_Roothan.o: solving_Roothan.cpp definitions.h
	c++ -c solving_Roothan.cpp

filling_matrices.o: filling_matrices.cpp definitions.h
	c++ -c filling_matrices.cpp

energy_orbitals.o: energy_orbitals.cpp definitions.h
	c++ -c energy_orbitals.cpp

Cs_evol.o: Cs_evol.cpp definitions.h
	c++ -c Cs_evol.cpp

der_matrix_elements.o: der_matrix_elements.cpp definitions.h
	c++ -c der_matrix_elements.cpp

nuclear_matrices.o: nuclear_matrices.cpp definitions.h
	c++ -c nuclear_matrices.cpp

X_evol.o: X_evol.cpp definitions.h
	c++ -c X_evol.cpp

CG_routine.o: CG_routine.cpp definitions.h
	c++ -c CG_routine.cpp

DFT_xc.o: DFT_xc.cpp definitions.h
	c++ -c DFT_xc.cpp

DFT_derivatives.o: DFT_derivatives.cpp definitions.h
	c++ -c DFT_derivatives.cpp

Adaptive.o: Adaptive.cpp definitions.h
	c++ -c Adaptive.cpp

file_handling.o: file_handling.cpp definitions.h
	c++ -c file_handling.cpp


main_H2.o: main_H2.cpp definitions.h
	c++ -c main_H2.cpp

clean :
	rm -f H2 main_H2.o matrix_elements.o solving_Roothan.o filling_matrices.o energy_orbitals.o Cs_evol.o der_matrix_elements.o nuclear_matrices.o X_evol.o CG_routine.o DFT_xc.o DFT_derivatives.o Adaptive.o file_handling.o