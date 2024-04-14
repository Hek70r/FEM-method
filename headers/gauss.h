#pragma once
#ifndef GAUSS_H
#define GAUSS_H

#include <vector>
using namespace std;

class gaussNode;
class gaussNode2D;

class gaussNode {
public:
	double x;
	double weight;

	gaussNode(double xk0, double weight0);
};

class gaussNode2D {
public:
	gaussNode x;
	gaussNode y;

	gaussNode2D(gaussNode x0, gaussNode y0);
};

class calcPoints {
public:
	int nodeNum;
	int pcNum;
	vector<gaussNode2D> points;

	calcPoints(int nodeNum0);
	void printPoints();
	gaussNode2D& operator[](int index);
};

extern vector<vector<gaussNode>> nodes;

double gaussIntegration1D(unsigned int nodeNum0, double(*f)(double));
double gaussIntegration2D(unsigned int k, double(*f)(double, double));

#endif // GAUSS_H


