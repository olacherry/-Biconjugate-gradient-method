#include <iostream>
#include <vector>
#include <cmath>

// умножение матрицы на вектор
std::vector<double> matVecMultiply(const std::vector<std::vector<double>>& A, const std::vector<double>& x) {
    int n = A.size();
    std::vector<double> result(n, 0.0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += A[i][j] * x[j];
        }
    }

    return result;
}

// вычитание векторов
std::vector<double> vecSubtract(const std::vector<double>& a, const std::vector<double>& b) {
    int n = a.size();
    std::vector<double> result(n, 0.0);

    for (int i = 0; i < n; i++) {
        result[i] = a[i] - b[i];
    }

    return result;
}

// 
double vecDotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    int n = a.size();
    double result = 0.0;

    for (int i = 0; i < n; i++) {
        result += a[i] * b[i];
    }

    return result;
}

// регение системы методом
std::vector<double> biconjugateGradient(const std::vector<std::vector<double>>& A, const std::vector<double>& b, int maxIterations, double tolerance) {
    int n = A.size();
    std::vector<double> x(n, 0.0);
    std::vector<double> r = vecSubtract(b, matVecMultiply(A, x));
    std::vector<double> r0 = r;
    std::vector<double> p = r;
    std::vector<double> p0 = p;

    for (int i = 0; i < maxIterations; i++) {
        std::vector<double> Ap = matVecMultiply(A, p);
        double alpha = vecDotProduct(r0, r) / vecDotProduct(Ap, p);
        std::vector<double> x_new = vecSubtract(x, vecScalarMultiply(alpha, p));
        std::vector<double> r_new = vecSubtract(r, vecScalarMultiply(alpha, Ap));

        if (vecDotProduct(r_new, r_new) < tolerance * tolerance) {
            return x_new;
        }

        double beta = vecDotProduct(r0, r_new) / vecDotProduct(r0, r0);
        std::vector<double> p_new = vecAdd(r_new, vecScalarMultiply(beta, vecSubtract(p, vecScalarMultiply(beta, Ap))));

        x = x_new;
        r = r_new;
        p = p_new;
    }

    return x;
}

int main() {
    // пример 
    std::vector<std::vector<double>> A = {{3.0, -1.0}, {-1.0, 2.0}};
    std::vector<double> b = {2.0, 1.0};

    std::vector<double> solution = biconjugateGradient(A, b, 1000, 1e-6);

    std::cout << "Solution: ";
    for (double x : solution) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    return 0;
}
