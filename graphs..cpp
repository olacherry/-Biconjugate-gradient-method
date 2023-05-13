#include <iostream>
#include <vector>
#include <cmath>
#include <matplotlibcpp.h>
#include <main.cpp>

namespace plt = matplotlibcpp;

// ...


std::vector<double> biconjugateGradient(const std::vector<std::vector<double>>& A, const std::vector<double>& b, int maxIterations, double tolerance) {
    // ...
    
    std::vector<double> residualNorms;
    
    for (int i = 0; i < maxIterations; i++) {
        // ...
        
        double residualNorm = std::sqrt(vecDotProduct(r_new, r_new));
        residualNorms.push_back(residualNorm);
        
        // ...
    }
    
    return residualNorms;
}

int main() {
    // ...

    std::vector<double> residualNorms = biconjugateGradient(A, b, 1000, 1e-6);

    // Plotting the convergence
    plt::plot(residualNorms);
    plt::xlabel("Iteration");
    plt::ylabel("Residual Norm");
    plt::show();

    return 0;
}
