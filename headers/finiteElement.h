#pragma once
#ifndef UNIVERSAL_ELEMENT_H
#define UNIVERSAL_ELEMENT_H
#include "global_data.h"
#include "matrix.h"
#include "Side.h"
#include <vector>

double N1(double ksi, double eta);
double N2(double ksi, double eta);
double N3(double ksi, double eta);
double N4(double ksi, double eta);
double dN1dETA(double ksi);
double dN2dETA(double ksi);
double dN3dETA(double ksi);
double dN4dETA(double ksi);
double dN1dKSI(double eta);
double dN2dKSI(double eta);
double dN3dKSI(double eta);
double dN4dKSI(double eta);
class Side;
class finiteElement;
class node;

class finiteElement {
public:
	
	int nodeNum;
	int pcNum;
	const int sfNum = 4;

	vector<node> coordinates;
	vector<vector<double>> Ni;
	vector<vector<double>> dNdKSI;
	vector<vector<double>> dNdETA;
	vector<vector<double>> dNdX;
	vector<vector<double>> dNdY;
	vector<double> detJ;
	vector<vector<vector<double>>> Hpc;
	vector<vector<double>> Hbc;
	vector<vector<double>> H;
	vector<Side> sides;
	vector<double> P;
	vector<vector<vector<double>>> Cpc;
	vector<vector<double>> C;

	finiteElement();
	finiteElement(globalData globalD, vector<node> coordinates0, int nodeNum0);
	void resizeFields();
	void fillNi();
	void fillEtaKsiDeratives(); // fills the matrixes dNdKSI and dNdETA
	void fillXYDerativesAndJacobians();
	void calculateHpc(double conductiwity);
	void calculateH();
	void calculateCpcAndC(double c, double ro);
	void fillSides(globalData gData);
	void calculateHbc();
	void calculateP();

	const void showdNdX();
	const void showdNdY();
	const void showdNdETA();
	const void showdNdKSI();
	const void showH();
	const void showC();
	const void showHbc();
	const void showP();
	const void showCoordinates();
};
#endif