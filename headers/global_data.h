#pragma once
#ifndef GLOBAL_DATA_H
#define GLOBAL_DATA_H

#include <vector>
#include <string>
#include <chrono>
using namespace std;
using namespace std::chrono;

class node;
class element;
class globalData;
class finiteElement;

class globalData {
public: 
	double SimulationTime;
	double SimulationStepTime;
	double Conductivity;
	double Alfa;
	double Tot;
	double InitialTemp;
	double Density;
	double SpecificHeat;
	unsigned int NodesNumber;
	unsigned int ElementsNumber;
	vector<node> nodes;
	vector<element> elements;

	vector<vector<double>> globalH;
	vector<vector<double>> globalC;
	vector<double> globalP;

	// NONSTATIONARY SOLUTION PARAMETERS
	vector<double> t0;
	vector<double> t1;
	vector<vector<double>> A;
	vector<double> b;
	vector<double> Solution;

	globalData();
	void read_data_from_file(string);
	void show();
	void aggregateGlobalH(vector<finiteElement> universalElements);
	void aggregateGlobalP(vector<finiteElement> universalElements);
	void aggregateGlobalC(vector<finiteElement> universalElements);

	void calculateA();
	void calculateB();
	void calculateSolution();
	void createParaViewFile(int step);
};

class node {
public:
	double x;
	double y;
	int BC;
	node(double x0 = 0, double y0 = 0, int BC0 = 0);
};

class element {
public:
	const int number_of_components = 4;
	vector<int> components;

	element();
};

#endif // GLOBAL_DATA_H