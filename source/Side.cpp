#include "Side.h"
#include "finiteElement.h"
#include "matrix.h"
#include "gauss.h"

Side::Side() {

}
Side::Side(globalData gData, node start0, node end0, int nodeNum0, int sideNum0) :
    start(start0), end(end0), nodeNum(nodeNum0), sideNum(sideNum0) {
    
    vector<vector<double>> newNi(nodeNum, vector<double>(4, 0.0));
    this->Ni = newNi;
    this->ksi.resize(nodeNum);
    this->eta.resize(nodeNum);
    this->weights.resize(nodeNum);

    vector<vector<double>> newHbcs(4, vector<double>(4, 0.0));
    this->Hbcs = newHbcs;
    vector<double> newPs(4, 0.0);
    this->Ps = newPs;

    double lenght = sqrt((end.x - start.x) * (end.x - start.x) + (end.y - start.y) * (end.y - start.y));
    this->detJ = lenght / 2.0;

    // Depending of which side it is, calculating calculus points coordinates
    int nodeIndex = nodeNum - 1;
    for (int i = 0; i < nodeNum; i++) {
        if (sideNum == 1) {
            this->ksi[i] = nodes[nodeIndex][i].x;
            this->eta[i] = -1;
            this->weights[i] = nodes[nodeIndex][i].weight;
        }
        else if (sideNum == 2) {
            this->ksi[i] = 1;
            this->eta[i] = nodes[nodeIndex][i].x;
            this->weights[i] = nodes[nodeIndex][i].weight;
        }
        else if (sideNum == 3) {
            this->ksi[i] = nodes[nodeIndex][i].x;
            this->eta[i] = 1;
            this->weights[i] = nodes[nodeIndex][i].weight;
        }
        else if (sideNum == 4) {
            this->ksi[i] = -1;
            this->eta[i] = nodes[nodeIndex][i].x;
            this->weights[i] = nodes[nodeIndex][i].weight;
        }
    }

    if (this->start.BC > 0 && this->end.BC > 0) {
        this->fillNi();
        this->calculateHbc(gData.Alfa);
        this->calculateP(gData.Alfa, gData.Tot);
    }
}
void Side::fillNi() {
    for (int i = 0; i < this->nodeNum; i++) {
        Ni[i][0] = N1(this->ksi[i], this->eta[i]);
        Ni[i][1] = N2(this->ksi[i], this->eta[i]);
        Ni[i][2] = N3(this->ksi[i], this->eta[i]);
        Ni[i][3] = N4(this->ksi[i], this->eta[i]);
    }
}
void Side::calculateHbc(double alfa) {
    vector<vector<double>> Hbcp;

    for (int i = 0; i < nodeNum; i++) {
        Hbcp = VectorVectorTMultiply(Ni[i]);
        MatrixMultiplyScalar(Hbcp, this->weights[i] * alfa);
        MatrixAdd(this->Hbcs, Hbcp);
    }
    MatrixMultiplyScalar(this->Hbcs, this->detJ);

}
void Side::calculateP(double alfa, double Tot) {

    vector<double> Pp;
    for (int i = 0; i < nodeNum; i++) {
        Pp = newVectorMultiplyScalar(Ni[i], Tot * weights[i]);
        VectorAdd(this->Ps, Pp);
    }
    VectorMulitplyScalar(this->Ps, alfa * detJ);
}