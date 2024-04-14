#pragma once
#ifndef SIDE_H
#define SIDE_H

#include <vector>
#include "global_data.h"
using namespace std;

class node;
class globalData;

class Side {
public:
	node start;
	node end;
	int nodeNum;
	int sideNum;
	double detJ;
	vector<double> eta;
	vector<double> ksi;
	vector<double> weights;
	vector<vector<double>> Ni;
	vector<vector<double>> Hbcs;
	vector<double> Ps;

	Side();
	Side(globalData gData, node start0, node end0, int nodeNum0, int sideNum0);
	void fillNi();
	void calculateHbc(double alfa);
	void calculateP(double alfa, double Tot);
};

#endif