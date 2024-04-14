#include "finiteElement.h"
#include "gauss.h"
#include <iostream>

finiteElement::finiteElement() {
    this->nodeNum = 0;
    this->pcNum = 0;
}
finiteElement::finiteElement(globalData globalD, vector<node> coordinates0, int nodeNum0 = 2)
: coordinates(coordinates0), nodeNum(nodeNum0)
{
    if (nodeNum < 1 || nodeNum > 5) {
        throw std::invalid_argument("Node number is not in [1, 5] range");
    }

    this->pcNum = nodeNum * nodeNum;

    this->resizeFields();

    this->fillNi();
    this->fillEtaKsiDeratives();
    this->fillXYDerativesAndJacobians();
    this->calculateHpc(globalD.Conductivity);
    this->calculateH();
    this->calculateCpcAndC(globalD.SpecificHeat, globalD.Density );
    this->fillSides(globalD);
    this->calculateHbc();
    this->calculateP();
}

void finiteElement::resizeFields() {
    // Initialization of Shape functions deratives with respect to KSI and ETA
    // Initialization of Shape functions deratives with respect to X and Y
    this->Ni.resize(pcNum);
    this->dNdKSI.resize(pcNum); 
    this->dNdETA.resize(pcNum);
    this->dNdX.resize(pcNum); 
    this->dNdY.resize(pcNum);
    for (int i = 0; i < this->pcNum; i++) {
        Ni[i].resize(sfNum);
        dNdKSI[i].resize(sfNum); 
        dNdETA[i].resize(sfNum);
        dNdX[i].resize(sfNum); 
        dNdY[i].resize(sfNum);
    }
    this->detJ.resize(pcNum);

    // Inicialization of array of Hp matrixes for every calculus point
    this->Hpc.resize(0);
    this->H.resize(sfNum);
    this->Cpc.resize(0);
    this->C.resize(sfNum);
    this->Hbc.resize(sfNum);
    this->P.resize(sfNum);
    for (int i = 0; i < sfNum; i++) {
        this->H[i].resize(sfNum);
        this->C[i].resize(sfNum);
        this->Hbc[i].resize(sfNum);
        this->P[i] = 0.0;
        for (int j = 0; j < sfNum; j++) {
            H[i][j] = 0.0;
            Hbc[i][j] = 0.0;
        }
    }
    
    // Inicialization of array of Sides (4) 
    this->sides.resize(sfNum);
    this->P.resize(sfNum);

}
void finiteElement::fillNi() {
    calcPoints points(this->nodeNum);
    for (int i = 0; i < this->pcNum; i++) {
        Ni[i][0] = N1(points[i].x.x, points[i].y.x);
        Ni[i][1] = N2(points[i].x.x, points[i].y.x);
        Ni[i][2] = N3(points[i].x.x, points[i].y.x);
        Ni[i][3] = N4(points[i].x.x, points[i].y.x);
    }
}

void finiteElement::fillEtaKsiDeratives() {
    calcPoints points(this->nodeNum);

    for (int i = 0; i < this->pcNum; i++) {
        dNdETA[i][0] = dN1dETA(points.points[i].x.x);
        dNdKSI[i][0] = dN1dKSI(points.points[i].y.x);
  
        dNdETA[i][1] = dN2dETA(points.points[i].x.x);
        dNdKSI[i][1] = dN2dKSI(points.points[i].y.x);

        dNdETA[i][2] = dN3dETA(points.points[i].x.x);
        dNdKSI[i][2] = dN3dKSI(points.points[i].y.x);

        dNdETA[i][3] = dN4dETA(points.points[i].x.x);
        dNdKSI[i][3] = dN4dKSI(points.points[i].y.x);
    }
}
void finiteElement::fillXYDerativesAndJacobians() {
    // Przejœcie po wszystkich pkt ca³kowania, bo dla ka¿dego jest osobny Jakobian, Macierz z dydETA -dydKSI -dxdETA dxdKSI
    for (int i = 0; i < this->pcNum; i++) {

        // Liczenie tej poszczególnych wartoœci tej macierzy  dydETA -dydKSI -dxdETA dxdKSI
        double dYdETA = 0.0, dYdKSI = 0.0, dXdETA = 0.0, dXdKSI = 0.0;
        for (int j = 0; j < sfNum; j++) { // sfNum = 4
            dYdETA += dNdETA[i][j] * this->coordinates[j].y;
            dYdKSI += dNdKSI[i][j] * this->coordinates[j].y;
            dXdETA += dNdETA[i][j] * this->coordinates[j].x;
            dXdKSI += dNdKSI[i][j] * this->coordinates[j].x;
        }

        double Jacobian;
        Jacobian = dXdKSI * dYdETA - dYdKSI * dXdETA;
        //cout << "Jacobian nr." << i << " = " << Jacobian << "\n";
        this->detJ[i] = Jacobian;

        for (int j = 0; j < sfNum; j++) { // sfNum = 4
            this->dNdX[i][j] = (1 / Jacobian) * (dYdETA * dNdKSI[i][j] + (-dYdKSI) * dNdETA[i][j]);
            this->dNdY[i][j] = (1 / Jacobian) * ((-dXdETA) * dNdKSI[i][j] + dXdKSI * dNdETA[i][j]);
        }
    }
}

void finiteElement::calculateHpc(double conductivity) {

    for (int i = 0; i < this->pcNum; i++) {
        vector<vector<double>> newHpc(this->sfNum, vector<double>(this->sfNum, 0.0)); // Zawsze 4x4
        vector<vector<double>> halfH = VectorVectorTMultiply(dNdX[i]);
        MatrixAdd(newHpc, halfH);
        halfH = VectorVectorTMultiply(dNdY[i]);
        MatrixAdd(newHpc, halfH);
        MatrixMultiplyScalar(newHpc, conductivity * detJ[i]);

        this->Hpc.push_back(newHpc);
    }
}
void finiteElement::calculateH() {
    calcPoints points(this->nodeNum);

    for (int i = 0; i < this->pcNum; i++) {
        double weight = points[i].x.weight * points[i].y.weight;
        MatrixMultiplyScalar(this->Hpc[i], weight);
        MatrixAdd(this->H, this->Hpc[i]);
    }
}

void finiteElement::calculateCpcAndC(double c, double ro) {
    vector<vector<double>> newCpc;
    for (int i = 0; i < pcNum; i++) {
        newCpc = VectorVectorTMultiply(Ni[i]);
        //printMatrix(partialC);
        MatrixMultiplyScalar(newCpc, c * ro * this->detJ[i]);
        this->Cpc.push_back(newCpc);
    }

    calcPoints points(nodeNum);
    for(int i = 0; i < pcNum; i++) {
        double weight = points[i].x.weight * points[i].y.weight;
        MatrixMultiplyScalar(this->Cpc[i], weight);
        MatrixAdd(this->C, this->Cpc[i]);
    }
}
void finiteElement::fillSides(globalData globalD) {
    for (int i = 0; i < sfNum; i++) {
        int nextIndex = (i + 1) % sfNum;
        sides[i] = Side(globalD, this->coordinates[i], this->coordinates[nextIndex], this->nodeNum, i + 1);
    }
}
void finiteElement::calculateHbc() {
    for (int i = 0; i < sides.size(); i++) {
        MatrixAdd(Hbc, sides[i].Hbcs);
    }
}
void finiteElement::calculateP() {
    for (int i = 0; i < sides.size(); i++) {
        VectorAdd(P, sides[i].Ps);
    }
}

    const void finiteElement::showdNdETA() {
        cout << "PRINTING UNIVERSAL ELEMENT, dNdETA\n";
        printMatrix(this->dNdETA);
    }
    const void finiteElement::showdNdKSI() {
        cout << "\nPRINTING UNIVERSAL ELEMENT, dNdKSI\n";
        printMatrix(this->dNdKSI);
    }
    const void finiteElement::showdNdX() {
        cout << "\nPRINTING UNIVERSAL ELEMENT, dNdX\n";
        printMatrix(this->dNdX);
    }
    const void finiteElement::showdNdY() {
        cout << "\nPRINTING UNIVERSAL ELEMENT, dNdX\n";
        printMatrix(this->dNdY);
    }
    const void finiteElement::showH() {
        printMatrix(this->H);
    }
    const void finiteElement::showC() {
        printMatrix(this->C);
    }
    const void finiteElement::showHbc() {
        printMatrix(this->Hbc);
    }
    const void finiteElement::showP() {
        printVector(this->P);
    }
    const void finiteElement::showCoordinates() {
        cout << "1 side [x, y] = [" << this->sides[0].start.x << ", " << this->sides[0].end.x << "]\n";
        cout << "2 side [x, y] = [" << this->sides[1].start.x << ", " << this->sides[1].end.x << "]\n";
        cout << "3 side [x, y] = [" << this->sides[2].start.x << ", " << this->sides[2].end.x << "]\n";
        cout << "4 side [x, y] = [" << this->sides[3].start.x << ", " << this->sides[3].end.x << "]\n";
    }

double N1(double ksi, double eta) {
    return (1.0 / 4.0) * (1 - ksi) * (1 - eta); 
}
double N2(double ksi, double eta) {
    return (1.0 / 4.0) * (1 + ksi) * (1 - eta);
}
double N3(double ksi, double eta) {
    return (1.0 / 4.0) * (1 + ksi) * (1 + eta);
}
double N4(double ksi, double eta) {
    return (1.0 / 4.0) * (1 - ksi) * (1 + eta);
}
double dN1dETA(double ksi) {
    return (-(1.0 / 4.0) * (1 - ksi));
}
double dN2dETA(double ksi) {
    return (-(1.0 / 4.0) * (1 + ksi));
}
double dN3dETA(double ksi) {
    return ((1.0 / 4.0) * (1 + ksi));
}
double dN4dETA(double ksi) {
    return ((1.0 / 4.0) * (1 - ksi));
}
double dN1dKSI(double eta) {
    return (-(1.0 / 4.0) * (1 - eta));
}
double dN2dKSI(double eta) {
    return ((1.0 / 4.0) * (1 - eta));
}
double dN3dKSI(double eta) {
    return ((1.0 / 4.0) * (1 + eta));
}
double dN4dKSI(double eta) {
    return (-(1.0 / 4.0) * (1 + eta));
}
