#include <iostream>
#include <vector>
#include <cmath>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;

// Function to solve a linear system using the biconjugate gradient method and return the number of iterations
int biconjugateGradient(const std::vector<std::vector<double>>& A, const std::vector<double>& b, double tolerance) {
    int n = A.size();
    std::vector<double> x(n, 0.0);
    std::vector<double> r = b;
    std::vector<double> r0 = r;
    std::vector<double> p = r;
    std::vector<double> p0 = p;
    
    int iterations = 0;

    while (true) {
        double alpha = 0.0;
        double beta = 0.0;

        std::vector<double> Ap(n, 0.0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Ap[i] += A[i][j] * p[j];
            }
        }

        double rr0 = 0.0;
        for (int i = 0; i < n; i++) {
            rr0 += r[i] * r0[i];
        }

        double pAp = 0.0;
        for (int i = 0; i < n; i++) {
            pAp += p[i] * Ap[i];
        }

        alpha = rr0 / pAp;

        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rr = 0.0;
        for (int i = 0; i < n; i++) {
            rr += r[i] * r[i];
        }

        double rr0new = 0.0;
        for (int i = 0; i < n; i++) {
            rr0new += r[i] * r0[i];
        }

        beta = rr0new / rr;

        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i];
        }

        iterations++;

        double residualNorm = std::sqrt(rr);
        if (residualNorm < tolerance) {
            break;
        }
    }

    return iterations;
}

int main() {
    // Generate data for different matrix sizes
    std::vector<int> matrixSizes = {10, 20, 30, 40, 50};
    std::vector<int> iterations;

    for (int size : matrixSizes) {
        std::vector<std::vector<double>> A(size, std::vector<double>(size, 0.0));
        std::vector<double> b(size, 1.0);
        
        // Generate a symmetric positive-definite matrix
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                A[i][j] = std::abs(i - j) + 1.0;
            }
        }

        int numIterations = biconjugateGradient(A, b, 1e-6);
        iterations.push_back(numIterations);
    }

    // Plotting the number of iterations versus matrix size
    plt::plot(matrixSizes, iterations);
    plt::xlabel("Matrix Size");
    plt::ylabel("Number of Iterations");
    plt::show();
    return 0;
}
