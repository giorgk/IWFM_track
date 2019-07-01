#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <vector>


template<typename T>
class iwfmMatrix2D {
public:
	iwfmMatrix2D();
	void readFromFile(std::string filename, bool zBase);
	void fromVector(std::vector<std::vector<T>> VEC, bool zBase);
	T get(int i, int j);

	int nrows();
	int ncols();


private:
	std::vector<std::vector<T>> MAT;
	bool zeroBase;
	int Nrow;
	int Ncol;
};

template<typename T>
iwfmMatrix2D<T>::iwfmMatrix2D(){}

template<typename T>
void iwfmMatrix2D<T>::readFromFile(std::string filename, bool zBase) {
	zeroBase = zBase;

	MAT.clear();

	std::ifstream datafile(filename.c_str());
	if (!datafile.good()) {
		std::cout << "Can't open the file" << filename << std::endl;
	}
	else {
		std::cout << "Reading " << filename << " ..." << std::endl;
		char buffer[1024];
		datafile.getline(buffer, 1024);
		std::istringstream inp(buffer);
		inp >> Nrow;
		inp >> Ncol;
		for (int i = 0; i < Nrow; ++i) {
			datafile.getline(buffer, 1024);
			std::istringstream inp(buffer);
			std::vector<T> tmp;
			T tmpd;
			for (int j = 0; j < Ncol; ++j) {
				inp >> tmpd;
				tmp.push_back(tmpd);
			}
			MAT.push_back(tmp);
		}
	}
	datafile.close();
}


template<typename T>
void iwfmMatrix2D<T>::fromVector(std::vector<std::vector<T>> VEC, bool zBase) {
	MAT = VEC;
	zeroBase = zBase;
	Nrow = static_cast<int>(MAT.size());
	Ncol = static_cast<int>(MAT[0].size());
}

template<typename T>
T iwfmMatrix2D<T>::get(int i, int j) {
	if (i < Nrow && j < Ncol)
		return MAT[i][j];
	else
		return static_cast<T>(-9999999);
}

template<typename T>
int iwfmMatrix2D<T>::nrows() {
	return Nrow;
}

template<typename T>
int iwfmMatrix2D<T>::ncols() {
	return Ncol;
}