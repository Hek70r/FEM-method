#include "global_data.h"
#include "matrix.h"
#include "finiteElement.h"
#include<iostream>
#include<iomanip>
#include<string>
#include<fstream>


// CONSTRUCTORS
globalData::globalData() {
    this->SimulationTime = 0;
    this->SimulationStepTime = 0;
    this->Conductivity = 0;
    this->Alfa = 0;
    this->Tot = 0;
    this->InitialTemp = 0;
    this->Density = 0;
    this->SpecificHeat = 0;
    this->NodesNumber = 0;
    this->ElementsNumber = 0;
    this->nodes.resize(0);
    this->elements.resize(0);
}
node::node(double x0, double y0, int BC0) {
    this->x = x0;
    this->y = y0;
    this->BC = BC0;
}
element::element() {
    this->components = vector<int>(this->number_of_components, 0.0);
}


// REDING FROM FILE
void globalData::read_data_from_file(string filePath) {

    fstream file(filePath);
    if (!file.is_open()) {
        throw runtime_error("Nie mo¿na otworzyæ pliku " + filePath);
    }

    string key;
    // Reading global data
    for (int i = 0; i < 10; i++) {
        if (file >> key) {
            if (key == "SimulationTime") {
                file >> SimulationTime;
            }
            else if (key == "SimulationStepTime") {
                file >> SimulationStepTime;
            }
            else if (key == "Conductivity") {
                file >> Conductivity;
            }
            else if (key == "Alfa") {
                file >> Alfa;
            }
            else if (key == "Tot") {
                file >> Tot;
            }
            else if (key == "InitialTemp") {
                file >> InitialTemp;
            }
            else if (key == "Density") {
                file >> Density;
            }
            else if (key == "SpecificHeat") {
                file >> SpecificHeat;
            }
            else if (key == "Nodes" && file >> key && key == "number") {
                file >> NodesNumber;
            }
            else if (key == "Elements" && file >> key && key == "number") {
                file >> ElementsNumber;
            }
        }
    }

    this->nodes     = vector<node>(this->NodesNumber, node());
    this->elements  = vector<element>(this->ElementsNumber, element());

    // reading nodes data
    string help_str;
    while (file >> key) {
        if (key == "*Node") {
            for (int i = 0; i < this->NodesNumber; i++) {
                file >> help_str; // number of node, does nothing
                file >> help_str; this->nodes[i].x = stod(help_str.substr(0, help_str.length() - 1)); // cutting of the comma (,) at the end and converting to double
                file >> help_str; this->nodes[i].y = stod(help_str); // there is no comma (,) in the last value in the line
            }
            break;
        }
    }

    // reading elements data
    while (file >> key) {
        if (key == "*Element," && file >> key) {
            for (int i = 0; i < this->ElementsNumber; i++) {
                file >> help_str; // number of element, does nothing
                file >> help_str; this->elements[i].components[0] = stoi(help_str.substr(0, help_str.length() - 1)); // cutting of the comma (,) at the end and converting to int
                file >> help_str; this->elements[i].components[1] = stoi(help_str.substr(0, help_str.length() - 1)); // cutting of the comma (,) at the end and converting to int
                file >> help_str; this->elements[i].components[2] = stoi(help_str.substr(0, help_str.length() - 1)); // cutting of the comma (,) at the end and converting to int
                file >> help_str; this->elements[i].components[3] = stoi(help_str); // there is no comma (,) in the last value in the line
            }
            break;
        }
    }
    vector<string> bc_vec(0);
    while (file >> key) {
        if (key == "*BC") {
            while (file >> key) {
                bc_vec.push_back(key);
            }
            // cutting of the comma (,) after every number except of the last one
            for (int i = 0; i < bc_vec.size() - 1; i++) {
                bc_vec[i] = bc_vec[i].substr(0, bc_vec[i].length() - 1);
            }
            //  adding BC to appropriate nodes
            for (int i = 0; i < bc_vec.size(); i++) {
                int nodeIndex = stoi(bc_vec[i]) - 1;
                this->nodes[nodeIndex].BC++;
            }
        }
    }
}
void globalData::show() {
    cout << "SimulationTime = " << this->SimulationTime << "\n";
    cout << "SimulationStepTime = " << this->SimulationStepTime << "\n";
    cout << "Conductivity = " << this->Conductivity << "\n";
    cout << "Alfa = " << this->Alfa << "\n";
    cout << "Tot = " << this->Tot << "\n";
    cout << "InitialTemp = " << this->InitialTemp << "\n";
    cout << "Density = " << this->Density << "\n";
    cout << "SpecificHeat = " << this->SpecificHeat << "\n";
    cout << "NodesNumber = " << this->NodesNumber << "\n";
    cout << "ElementsNumber = " << this->ElementsNumber << "\n";

    cout << "nodes:\n";
    for (int i = 0; i < this->NodesNumber; i++) {
        cout << i+1 << ".\t[x, y, BC] = " << "[" << this->nodes[i].x << ", " << this->nodes[i].y << ", " << this->nodes[i].BC << "]\n";
    }

    cout << "elements:\n";
    for (int i = 0; i < this->ElementsNumber; i++) {
        cout << i+1 << ".\t[id1, id2, id3, id4] = " << "[" << this->elements[i].components[0] << ", " << this->elements[i].components[1] << ", " << 
            this->elements[i].components[2] << ", " <<  this->elements[i].components[3] << "]\n";
    }
}

// AGGREGATION
void globalData::aggregateGlobalH(vector<finiteElement> universalElements) {

    globalH = vector<vector<double>>(NodesNumber, vector<double>(NodesNumber, 0.0));

    for (int i = 0; i < ElementsNumber; i++) {
        for (int j = 0; j < 4; j++) {
            int index1 = this->elements[i].components[j] - 1; //indexowanie od 0, a w pliku od 1
            for (int k = 0; k < 4; k++) {
                int index2 = this->elements[i].components[k] - 1;
                double val1 = universalElements[i].H[j][k];
                double val2 = universalElements[i].Hbc[j][k];

                this->globalH[index1][index2] += (val1 + val2);
            }
        }
    }

    //cout << "\nPO OBLICZENIU, GLOBAL H:\n";
    //printMatrix(globalH);
}
void globalData::aggregateGlobalP(vector<finiteElement> universalElements) {
    globalP = vector<double>(NodesNumber, 0.0);
    
    for (int i = 0; i < ElementsNumber; i++) {
        for (int j = 0; j < 4; j++) {
            int index = this->elements[i].components[j] - 1;
            this->globalP[index] += universalElements[i].P[j];
        }
    }

    //printVector(this->globalP);
}
void globalData::aggregateGlobalC(vector<finiteElement> universalElements) {

    globalC = vector<vector<double>>(NodesNumber, vector<double>(NodesNumber, 0.0));

    for (int i = 0; i < ElementsNumber; i++) {
        for (int j = 0; j < 4; j++) {
            int index1 = this->elements[i].components[j] - 1;
            for (int k = 0; k < 4; k++) {
                int index2 = this->elements[i].components[k] - 1;
                double val1 = universalElements[i].C[j][k];

                this->globalC[index1][index2] += val1;
            }
        }
    }

    //cout << "\nPO OBLICZENIU, GLOBAL C:\n";
    //printMatrix(globalC);
}

// SOLUTION
void globalData::calculateA() {
    this->A = vector<vector<double >>(NodesNumber, vector<double>(NodesNumber, 0.0));
    MatrixAdd(A, this->globalH);
    MatrixAdd(A, newMatrixMultiplyScalar(this->globalC, 1.0 / SimulationStepTime));
    //cout << "\n [H] + [C]/dTau:";
    //printMatrix(A);
}
void globalData::calculateB() {
    this->b = vector<double>(NodesNumber, 0.0);
    for (int j = 0; j < NodesNumber; j++) {
        b[j] = globalP[j];
        for (int k = 0; k < NodesNumber; k++) {
            b[j] += (globalC[j][k] / SimulationStepTime) * t0[k];
        }
    }
    //cout << "\n -([C]/dTau){t0} + {P}:";
    //printVector(b);
}
void globalData::calculateSolution() {
    this->t0 = vector<double>(NodesNumber, InitialTemp);
    this->t1 = vector<double>(NodesNumber, 0.0);
    this->Solution = vector<double>(NodesNumber, 0.0);
    int simulationSteps = SimulationTime / SimulationStepTime;

    for (int i = 0; i < simulationSteps; i++) {
        this->calculateA();
        this->calculateB();


        auto start = high_resolution_clock::now();
        this->Solution = SolveLinearSystem(this->A, this->b);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Czas liczenia ROZWIAZANIA: " << duration.count() / 1e+6 << " sekund" << endl;

        double max = findMaxInVector(this->Solution);
        double min = findMinInVector(this->Solution);

        cout << "dTau = " << SimulationStepTime * (i + 1) 
            << "\tt_min = " << std::setprecision(10) << min 
            << "\tt_max = " << std::setprecision(10) << max << "\n";

        t0 = this->Solution;
        t1 = vector<double>(NodesNumber, 0.0);
        createParaViewFile(i + 1);
    }
}

void globalData::createParaViewFile(int step) {
    string filename = "./ParaView/Foo" + to_string(step) + ".vtk";
    fstream file(filename, std::ios::out | std::ios::trunc);

    if (!file.is_open()) {
        cerr << "B³¹d otwarcia pliku: " << filename << endl;
        return;
    }

    file << "# vtk DataFile Version 2.0\n";
    file << "Unstructured Grid Example\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    file << "\nPOINTS " << NodesNumber << " float\n";
    for (int i = 0; i < NodesNumber; i++) {
        file << this->nodes[i].x << " " << this->nodes[i].y << " 0\n";
    }

    file << "\nCELLS " << ElementsNumber << " " << ElementsNumber * 5 << "\n";
    for (int i = 0; i < ElementsNumber; i++) {
        file << "4 " << elements[i].components[0] - 1 << " ";
        file << elements[i].components[1] - 1 << " ";
        file << elements[i].components[2] - 1 << " ";
        file << elements[i].components[3] - 1 << "\n";
    }

    file << "\nCELL_TYPES " << ElementsNumber << "\n";
    for (int i = 0; i < ElementsNumber; i++) {
        file << "9\n";
    }

    file << "\nPOINT_DATA " << NodesNumber << "\n";
    file << "SCALARS Temp float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < NodesNumber; i++) {
        file << this->Solution[i] << "\n";
    }

    file.close();
}
