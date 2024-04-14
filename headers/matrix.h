#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include <vector>\
using namespace std;

void MatrixAdd(std::vector <std::vector<double>>& a, std::vector <std::vector<double>> b);
void VectorAdd(std::vector<double>& a, std::vector<double> b);

void MatrixMultiplyScalar(std::vector < std::vector<double>>& a, double x);
void VectorMulitplyScalar(std::vector<double>& a, double x);

std::vector<double> newVectorMultiplyScalar(std::vector<double> a, double x);
std::vector<std::vector<double>> newMatrixMultiplyScalar(std::vector<std::vector<double>> a, double x);
std::vector <std::vector<double>>VectorVectorTMultiply(std::vector<double> a);
std::vector<double> MatrixVectorMultiply(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector);

void printMatrix(std::vector<std::vector<double>> a);
void printVector(std::vector<double> a);

vector<double> solveLinearEquation(const vector<vector<double>>& A, const vector<double>& b);
std::vector<double> SolveLinearSystem(const std::vector<std::vector<double>> A, const std::vector<double> b);
std::vector<std::vector<double>> InvertMatrix2x2(const std::vector<std::vector<double>>& matrix);

double findMaxInVector(vector<double> vec);
double findMinInVector(vector<double> vec);

#endif
