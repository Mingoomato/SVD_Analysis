#include <iostream>
#include <clocale>
#include <locale>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
using namespace std;

int SVD(int rows, int cols, const vector<vector<int>>& matrix) {
    // 1. Transpose 구하기 (전치 행렬 값 대입)
    vector<vector<int>> transposed_matrix(cols, vector<int>(rows, 0));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            transposed_matrix[j][i] = matrix[i][j];
        }
    }

    // 2. A^T * A 곱셈 수행 (결과 행렬 U의 크기는 cols x cols)
    vector<vector<int>> U(cols, vector<int>(cols, 0));
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < cols; j++) {
            for (int k = 0; k < rows; k++) {
                U[i][j] += transposed_matrix[i][k] * matrix[k][j];
            }
        }
    }

    // 3. Eigen value of U
    vector<vector<int>> I(cols, vector<int>(cols, 0));
    double lambda;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (i == j) {
                I[i][j] = 1;
                U[i][j] -= lambda * I[i][j];
            }
            else
                I[i][j] = 0;
        }
    }

    // 4. Determinant via Gaussian Elimination of U.
    int pivot = 0;
	int max_pivot = 0;
	int determinant = 1;
	int sign = 1; // Track the number of row swaps for determinant sign adjustment
    for (int i = 0; i < cols; i++) {
		pivot = U[i][i];
        for (int j = i + 1; j < cols; j++) {
            if (pivot > U[j][i])
                continue;
            else if (pivot < U[j][i]) {
                pivot = U[j][i];
                swap(U[i], U[j]);
				sign *= -1; // Row swap changes the sign of the determinant
            }
        } // ---> Pivot 순으로 정렬.
		max_pivot = U[i][i];
        for (int k = i + 1; k < cols; k++) {
			double factor = static_cast<double>(U[k][i]) / max_pivot;
			for (int j = i; j < cols; j++) {
				U[k][j] -= static_cast<double>(factor * U[i][j]);
			}
        }
        
		determinant *= U[i][i];
    }
    
	determinant *= sign; // Adjust determinant sign based on row swaps

    return determinant;
}

int main() {
    // ensure wide/locale-aware output on Windows consoles
    setlocale(LC_ALL, "");
    // 1. 2D vector to store data (1218 rows, 7 columns)
    const int rows = 1218;
    const int cols = 7;
    vector<vector<int>> matrix(rows, vector<int>(cols));

    // 2. Open CSV file (use ASCII filename)
    // Rename or copy your CSV to "lotto_history.csv" or change the filename below.
    ifstream file("C:/Users/jmkjm/source/repos/Lottery_Analysis/Lottery_Analysis/lottowinnumber.csv");
    if (!file.is_open()) {
        cerr << "Failed to open the file. Check the path and file name!" << endl;
        return 1;
    }

    string line;
    // Skip header line
    getline(file, line);

    int currentRow = 0;
    while (getline(file, line) && currentRow < rows) {
        stringstream ss(line);
        string cell;

        // CSV parsing by comma
        // Columns C..I correspond to indices 2..8
        for (int i = 0; i < 9; ++i) {
            getline(ss, cell, ',');
            if (i >= 2 && i <= 8) {
                try {
                    matrix[currentRow][i - 2] = stoi(cell);
                }
                catch (...) {
                    matrix[currentRow][i - 2] = 0; // non-numeric -> 0
                }
            }
        }
        currentRow++;
    }

	cout << "determinant of Lottery Data " << SVD(rows, cols, matrix);

    return 0;
}