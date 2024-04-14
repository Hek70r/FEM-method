#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

void MatrixAdd(std::vector <std::vector<double>>& a, std::vector <std::vector<double>> b) {
    const int n = a.size();     // a rows
    const int m = a[0].size();  // a cols
    const int p = b.size();  // b rows
    const int k = b[0].size();  // b cols

    if (n != p || m != k) {
        std::cout << "Not valid matrix sizes\n";
    }

    for (auto i = 0; i < n; i++) {
        for (auto j = 0; j < m; j++) {
            a[i][j] += b[i][j];
        }
    }

}

void VectorAdd(std::vector<double>& a, std::vector<double> b) {
    for (int i = 0; i < a.size(); i++) {
        a[i] += b[i];
    }
}

void MatrixMultiplyScalar(std::vector < std::vector<double>>& a, double x) {
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < a[i].size(); j++) {
            a[i][j] *= x;
        }
    }
}

void VectorMulitplyScalar(std::vector<double>& a, double x) {
    for (int i = 0; i < a.size(); i++) {
        a[i] *= x;
    }
}

void printMatrix(std::vector<std::vector<double>> a) {
    std::cout << "{";
    for (int i = 0; i < a.size(); i++) {
        std::cout << "\n\t[";
        for (int j = 0; j < a[i].size(); j++) {
            std::cout << std::setprecision(6) << a[i][j] << ",   ";
        }
        std::cout << "]\n";
    }
    std::cout << "}";
}

void printVector(std::vector<double> a) {
    std::cout << "\n{\t";
    for (int i = 0; i < a.size(); i++) {
        std::cout << a[i] << "\t";
    }
    std::cout << "}\n";
}


std::vector <std::vector<double>> VectorVectorTMultiply(std::vector<double> a) {
    int size = a.size();
    std::vector<std::vector<double>> result(size, std::vector<double>(size, 0.0));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result[i][j] = a[i] * a[j];
        }
    }
    return result;
}

std::vector<double> MatrixVectorMultiply(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    if (cols != vector.size()) {
        std::cout << "Not valid matrix or vector sizes\n";
    }

    std::vector<double> result(rows, 0.0);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

std::vector<double> newVectorMultiplyScalar(std::vector<double> a, double x) {
    std::vector<double> result(a.size(), 0.0);
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] * x;
    }
    return result;
}

std::vector<std::vector<double>> newMatrixMultiplyScalar(std::vector<std::vector<double>> a, double x) {
    std::vector<std::vector<double>> result(a.size(), std::vector<double>(a.size(), 0.0));
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < a[i].size(); j++) {
            result[i][j] = a[i][j] * x;
        }
    }
    return result;
}

vector<double> solveLinearEquation(const vector<vector<double>>& A, const vector<double>& b) {
    vector<vector<double>> augmentedMatrix = A;
    for (size_t i = 0; i < A.size(); i++) {
        augmentedMatrix[i].push_back(b[i]);
    }

    size_t n = augmentedMatrix.size();
    for (size_t i = 0; i < n; i++) {
        double pivot = augmentedMatrix[i][i];
        for (size_t j = i; j < n + 1; j++) {
            augmentedMatrix[i][j] /= pivot;
        }

        for (size_t k = 0; k < n; k++) {
            if (k != i) {
                double factor = augmentedMatrix[k][i];
                for (size_t j = i; j < n + 1; j++) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }
    }

    vector<double> solution;
    for (size_t i = 0; i < n; i++) {
        solution.push_back(augmentedMatrix[i][n]);
    }

    return solution;
}

// Funkcja do rozwi¹zywania uk³adu równañ Ax + b = 0
std::vector<double> SolveLinearSystem(std::vector<std::vector<double>> A, std::vector<double> b) {
    int size = A.size();

    // Rozszerz macierz A o wektor b
    for (int i = 0; i < size; i++) {
        A[i].push_back(b[i]);
    }

    // Eliminacja wspó³czynników poni¿ej diagonalii
    for (int i = 0; i < size - 1; i++) {
        for (int k = i + 1; k < size; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j <= size; j++) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }

    // Ustalanie wartoœci nieznanych (wsteczna substitucja)
    std::vector<double> x(size, 0.0);
    for (int i = size - 1; i >= 0; i--) {
        x[i] = A[i][size];
        for (int j = i + 1; j < size; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    return x;
}


double findMaxInVector(vector<double> vec) {
    double max = vec[0];
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i] > max) {
            max = vec[i];
        }
    }
    return max;
}
double findMinInVector(vector<double> vec) {
    double min = vec[0];
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i] < min) {
            min = vec[i];
        }
    }
    return min;
}