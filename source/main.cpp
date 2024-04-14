#include "global_data.h"
#include "finiteElement.h"
#include "gauss.h"
#include "matrix.h"
#include "NonstationarySolution.h"
#include <iomanip>
#include <iostream>
#include <vector>5

const string INPUT_DATA_PATH = "./data/Test1_4_4.txt";

int main() {
	try {
		globalData Data;
		Data.read_data_from_file(INPUT_DATA_PATH);
		//Data.show();

		int CalcNodesNumber = 4;

		std::vector<finiteElement> Elements(0);
		for (int i = 0; i < Data.elements.size(); i++) {

			vector<node> currentNodes(4, 0.0);
			for (int j = 0; j < 4; j++) {
				int nodeIndex = Data.elements[i].components[j] - 1;
				currentNodes[j].x = Data.nodes[nodeIndex].x;
				currentNodes[j].y = Data.nodes[nodeIndex].y;
				currentNodes[j].BC = Data.nodes[nodeIndex].BC;
			}
			finiteElement newElement(Data, currentNodes, CalcNodesNumber);
			Elements.push_back(newElement);
		}

		Data.aggregateGlobalH(Elements);
		Data.aggregateGlobalP(Elements);
		Data.aggregateGlobalC(Elements);
		Data.calculateSolution();
	}
	catch (const std::exception& e) {
		cerr << "\nEXCEPTION!\n" << e.what();
	}

	return 0;
}