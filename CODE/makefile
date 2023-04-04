HF_H2: matrix_elements.o solving_Roothan.o filling_matrices.o energy_orbitals.o main_HF_H2.o
	c++ -o HF_H2 matrix_elements.o solving_Roothan.o filling_matrices.o energy_orbitals.o main_HF_H2.o -lm -lgsl -lgslcblas

matrix_elements.o: matrix_elements.cpp definitions.h
	c++ -c matrix_elements.cpp

solving_Roothan.o: solving_Roothan.cpp definitions.h
	c++ -c solving_Roothan.cpp

filling_matrices.o: filling_matrices.cpp definitions.h
	c++ -c filling_matrices.cpp

energy_orbitals.o: energy_orbitals.cpp definitions.h
	c++ -c energy_orbitals.cpp

main_HF_H2.o: main_HF_H2.cpp definitions.h
	c++ -c main_HF_H2.cpp

clean :
	rm -f HF_H2 main_HF_H2.o matrix_elements.o solving_Roothan.o filling_matrices.o energy_orbitals.o