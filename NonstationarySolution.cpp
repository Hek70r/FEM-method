#include "NonstationarySolution.h"
#include "global_data.h"
#include "matrix.h"
#include <iostream>
#include <fstream>
#include <string>

NonstationarySolution::NonstationarySolution(globalData gData) {
	this->NodesNumber = gData.NodesNumber;
	this->SimulationStepTime = gData.SimulationStepTime;
	this->SimulationTime = gData.SimulationTime;
	this->t0 = vector<double>(NodesNumber, gData.InitialTemp);
	this->t1 = vector<double>(NodesNumber, 0.0);
	this->H = gData.globalH;
	this->P = gData.globalP;
	this->C = gData.globalC;

	//this->A = vector<vector<double >>(NodesNumber, vector<double>(NodesNumber, 0.0));
	//this->b = vector<double>(NodesNumber, 0.0);
	//this->Solution = vector<double>(NodesNumber, 0.0);
}

void NonstationarySolution::calculateA() {
	this->A = vector<vector<double >>(NodesNumber, vector<double>(NodesNumber, 0.0));
	MatrixAdd(A, this->H);
	MatrixAdd(A, newMatrixMultiplyScalar(this->C, 1.0 / SimulationStepTime));
	//cout << "\n [H] + [C]/dTau:";
	//printMatrix(A);
}

void NonstationarySolution::calculateB() {
	//VectorAdd(b, P);
	//vector<double> help1 = MatrixVectorMultiply(newMatrixMultiplyScalar(this->C, 1.0 / dTau), this->t0);
	//VectorMulitplyScalar(help1, -1);
	//VectorAdd(b, help1);

	this->b = vector<double>(NodesNumber, 0.0);
	for (int j = 0; j < NodesNumber; j++) {
		b[j] = P[j];
		for (int k = 0; k < NodesNumber; k++) {
			b[j] += (C[j][k] / SimulationStepTime) * t0[k];
		}
	}
	//cout << "\n -([C]/dTau){t0} + {P}:";
	//printVector(b);
}

void NonstationarySolution::calculateSolution() {
	this->Solution = vector<double>(NodesNumber, 0.0);
	int simulationSteps = SimulationTime / SimulationStepTime;

	for (int i = 0; i < simulationSteps; i++) {
		this->calculateA();
		this->calculateB();

		this->Solution = SolveLinearSystem(this->A, this->b);
		double max = findMaxInVector(this->Solution);
		double min = findMinInVector(this->Solution);

		cout << "dTau = " << SimulationStepTime * (i+1) << "\tt_min = " << min << "\tt_max = " << max << "\n";

		t0 = this->Solution;
		t1 = vector<double>(NodesNumber, 0.0);

		createParaViewFile(i+1);
	}
	//cout << "\n SOLUTION:";
	//printVector(this->Solution);
}

void NonstationarySolution::createParaViewFile(int step) {
	string filename = "./ParaView/Foo.txt" + to_string(step);
	fstream file(filename, std::ios::out | std::ios::trunc);

	if (!file.is_open()) {
		cerr << "B³¹d otwarcia pliku: " << filename << endl;
		return;
	}

	file << "# vtk DataFile Version 2.0\n";
	file << "Unstructured Grid Example\n";
	file << "ASCII\n";
	file << "DATASET UNSTRUCTURED_GRID\n\n";

	file << "Points " << NodesNumber << " float\n";
	//for(

	file.close();
}