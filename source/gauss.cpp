#include "gauss.h"
#include <vector>
#include <iostream>

// Nodes classes =======================================================

gaussNode::gaussNode(double xk0 = 0, double weight0 = 0) {
    this->x = xk0;
    this->weight = weight0;
};


gaussNode2D::gaussNode2D(gaussNode x0 = gaussNode(), gaussNode y0 = gaussNode()) {
	this->x = x0;
	this->y = y0;
}
//======================================================================


// Certain nodes tables ===============================================
vector<vector<gaussNode>> nodes{
	{ 
		gaussNode(0.0, 2.0)																						// for 1 node
	},
	{ 
		gaussNode(-1.0 / sqrt(3.0), 1.0),																			// for 2 nodes
		gaussNode(+1.0 / sqrt(3.0), 1.0)
	},												
	{ 
		gaussNode(-sqrt(3.0 / 5.0), 5.0 / 9.0),																	// for 3 nodes
		gaussNode(0.0, 8.0 / 9.0),
		gaussNode(sqrt(3.0 / 5.0), 5.0 / 9.0) 
	},
	{ 
		gaussNode(-sqrt(((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0))), (18.0 - sqrt(30.0)) / 36.0),				// for 4 nodes																					// for 4 nodes						
		gaussNode(-sqrt(((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0))), (18.0 + sqrt(30.0)) / 36.0),
		gaussNode(+sqrt(((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0))), (18.0 + sqrt(30.0)) / 36.0),
		gaussNode(+sqrt(((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0))), (18.0 - sqrt(30.0)) / 36.0)
	},			
	{ 
		gaussNode(-((1.0 / 3.0) * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0))), (322.0 - 13.0 * sqrt(70.0)) / 900.0),		// for 5 nodes
		gaussNode(-((1.0 / 3.0) * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0))), (322.0 + 13.0 * sqrt(70.0)) / 900.0),
		gaussNode(0.0, 128.0/225.0),
		gaussNode(+((1.0 / 3.0) * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0))), (322.0 + 13.0 * sqrt(70.0)) / 900.0),
		gaussNode(+((1.0 / 3.0) * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0))), (322.0 - 13.0 * sqrt(70.0)) / 900.0),
	}
};
//======================================================================

// [] operator =========================================================
gaussNode2D& calcPoints::operator[](int index) {
	if (index >= 0 && index < this->points.size()) {
		return points[index];
	}
	else {
		throw std::out_of_range("Index out of range");
	}
}
//======================================================================

// Construcotr initializing the 2D points made of gauss nodes ==========
calcPoints::calcPoints(int nodeNum0 = 2) {
	this->nodeNum = nodeNum0;
	this->pcNum = nodeNum0 * nodeNum0;
	int nodeIndex = nodeNum0 - 1;

	// Space allocation
	this->points.resize(this->pcNum);
	for (gaussNode2D el : this->points) {
		el = gaussNode2D();
	}
	
	// creating points in the following order:
	//	7	8	9
	//	4	5	6
	//	1	2	3 etc.
	int whichXvalue = 0; // na planszy wartoœci na osi x to powtarzaj¹ce siê node'y
	int whichYvalue = 0; // na planszy wartoœci na osi y to 
	for (int i = 0; i < this->pcNum; i++) {
		this->points[i].x = nodes[nodeIndex][whichXvalue];
		this->points[i].y = nodes[nodeIndex][whichYvalue];

		// whichXvalue dodawane co iteracje, a zerowane co iloœæ nodeów
		whichXvalue++;
		if (whichXvalue >= nodeNum0) { 
			whichXvalue = 0; 
		}
		// whichTvlue jest iterowane co iloœæ nodeów
		if ((i + 1) % nodeNum0 == 0) {
			whichYvalue++;
		}
	}
}
//======================================================================


// 2D points printing method ===========================================
void calcPoints::printPoints() {
	for (int i = 0; i < this->pcNum; i++) {
		cout << "[" << this->points[i].x.x << ", " << this->points[i].y.x << "]\t";
		if ((i + 1) % this->nodeNum == 0)
			cout << "\n";
	}
}
//======================================================================


double gaussIntegration1D(unsigned int k, double(*f)(double)) {
	if (k < 1 || k > 5) {
		throw std::invalid_argument("Node number is not in [1, 5] range");
	}
	double result = 0.0;
	int nodeIndex = k - 1;
	for (int i = 0; i < nodes[nodeIndex].size(); i++) {
		result += nodes[nodeIndex][i].weight * f(nodes[nodeIndex][i].x);
	}
	return result;
}

double gaussIntegration2D(unsigned int k, double(*f)(double, double)) {
	if (k < 1 || k > 5) {
		throw std::invalid_argument("Node number is not in [1, 5] range");
	}
	double result = 0.0;
	int nodeIndex = k - 1;

	for (int i = 0; i < nodes[nodeIndex].size(); i++) {
		for (int j = 0; j < nodes[nodeIndex].size(); j++) {
			double x = nodes[nodeIndex][i].x;
			double y = nodes[nodeIndex][j].x;
			double x_weight = nodes[nodeIndex][i].weight;
			double y_weight = nodes[nodeIndex][j].weight;

			result += x_weight * y_weight * f(x, y);
		}
	}

	return result;
}