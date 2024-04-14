#pragma once
#ifndef NONSTATIONARY_SOLUTION_H
#define NONSTATIONARY_SOLUTION_H

#include <vector>
using namespace std;

class globalData;

class NonstationarySolution {
public:
	double SimulationStepTime;
	double SimulationTime;
	double NodesNumber;
	vector<double> t0;
	vector<double> t1;
	vector<double> P;
	vector<vector<double>> H;
	vector<vector<double>> C;
	// Ax * b = 0
	// A = H + C/dTau
	// x = t1
	// b = -(C/dTau) * t0 * P
	vector<vector<double>> A;
	vector<double> b;
	vector<double>Solution;

	NonstationarySolution(globalData gData);
	void calculateA();
	void calculateB();
	void calculateSolution();
	void createParaViewFile(int step);
};

#endif