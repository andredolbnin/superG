#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream> //
#include <stdlib.h>
#include <string>

using namespace std;


int M; // по x
int N; // по y

double h1;
double h2;


double u(double x, double y) {
    return sqrt(4.0 + x * y);
}

double q(double x, double y) {
    return x + y;
}

double k(double x, double y) {
    return 4.0 + x + y;
}

double F(double x, double y) {
    return - (1.0 / (4.0 * (4.0 + x * y) * sqrt(4.0 + x * y))) * (y * ((x - 4.0) * y - y * y + 8.0) + x * ((y - 4.0) * x - x * x + 8.0)) + q(x, y) * u(x, y);
}

double psi_L(double y) {
    return - (4.0 + y) * y / 4.0 + 2.0;
}

double psi_R(double y) {
    return (8.0 + y) * y / (4.0 * sqrt(1.0 + y)) + 2.0 * sqrt(1.0 + y);
}

double psi_B(double x) {
    return - (4.0 + x) * x / 4.0;
}

double psi_T(double x) {
    return (7.0 + x) * x / (2.0 * sqrt(4.0 + 3.0 * x));
}

double w_x(double* w, int i, int j) {
    return (w[(i + 1) * (N + 1) + j] - w[i * (N + 1) + j]) / h1;
}

double w_x_back(double* w, int i, int j) {
    return (w[i * (N + 1) + j] - w[(i - 1) * (N + 1) + j]) / h1;
}

double w_y(double* w, int i, int j) {
    return (w[i * (N + 1) + j + 1] - w[i * (N + 1) + j]) / h2;
}

double w_y_back(double* w, int i, int j) {
    return (w[i * (N + 1) + j] - w[i * (N + 1) + j - 1]) / h2;
}

void fill(double* A, double* w, double * B, bool fillRight = false)
{
    for (int i = 0; i <= M; ++i) {
        for (int j = 0; j <= N; ++j) {
            if (i == 0 && j == 0) { // (A1, B1)
                double a = k((1.0 - 0.5) * h1, j * h2) * w_x_back(w, 1, j);
                double b = k(i * h1, (1.0 - 0.5) * h2) * w_y_back(w, i, 1);
                A[i * (N + 1) + j] = -(2 / h1) * a - (2 / h2) * b + (q(i * h1, j * h2) + 2 / h1) * w[i * (N + 1) + j];
                if (fillRight) B[i * (N + 1) + j] = F(i * h1, j * h2) + (2 / h1 + 2 / h2) * (psi_L(j * h2) + psi_B(i * h1)) / 2;

            } else if (i == M && j == 0) { // (A2, B1)
                double a = k((i - 0.5) * h1, j * h2) * w_x_back(w, i, j);
                double b = k(i * h1, (1.0 - 0.5) * h2) * w_y_back(w, i, 1);
                A[i * (N + 1) + j] = (2 / h1) * a - (2 / h2) * b + (q(i * h1, j * h2) + 2 / h1) * w[i * (N + 1) + j];
                if (fillRight) B[i * (N + 1) + j] = F(i * h1, j * h2) + (2 / h1 + 2 / h2) * (psi_B(i * h1) + psi_R(j * h2)) / 2;

            } else if (i == M && j == N) { // (A2, B2)
                double a = k((i - 0.5) * h1, j * h2) * w_x_back(w, i, j);
                double b = k(i * h1, (j - 0.5) * h2) * w_y_back(w, i, j);
                A[i * (N + 1) + j] = (2 / h1) * a + (2 / h2) * b + (q(i * h1, j * h2) + 2 / h1) * w[i * (N + 1) + j];
                if (fillRight) B[i * (N + 1) + j] = F(i * h1, j * h2) + (2 / h1 + 2 / h2) * (psi_R(j * h2) + psi_T(i * h1)) / 2;

            } else if (i == 0 && j == N) { // (A1, B2)
                double a = k((1.0 - 0.5) * h1, j * h2) * w_x_back(w, 1.0, j);
                double b = k(i * h1, (j - 0.5) * h2) * w_y_back(w, i, j);
                A[i * (N + 1) + j] = -(2 / h1) * a + (2 / h2) * b + (q(i * h1, j * h2) + 2 / h1) * w[i * (N + 1) + j];
                if (fillRight) B[i * (N + 1) + j] = F(i * h1, j * h2) + (2 / h1 + 2 / h2) * (psi_T(i * h1) + psi_L(j * h2)) / 2;

            } else if (i == 0) { // левая граница
                double a = k((1.0 - 0.5) * h1, j * h2) * w_x_back(w, 1, j);
                double b = k(i * h1, (j + 0.5) * h2) * w_y(w, i, j) - k(i * h1, (j - 0.5) * h2) * w_y_back(w, i, j);
                A[i * (N + 1) + j] = -(2 / h1) * a + (q(i * h1, j * h2) + 2 / h1) * w[i * (N + 1) + j] - b / h2;
                if (fillRight) B[i * (N + 1) + j] = F(i * h1, j * h2) + (2 / h1) * psi_L(j * h2);

            } else if (i == M) { // правая граница
                double a = k((i - 0.5) * h1, j * h2) * w_x_back(w, i, j);
                double b = k(i * h1, (j + 0.5) * h2) * w_y(w, i, j) - k(i * h1, (j - 0.5) * h2) * w_y_back(w, i, j);
                A[i * (N + 1) + j] = (2 / h1) * a + (q(i * h1, j * h2) + 2 / h1) * w[i * (N + 1) + j] - b / h2;
                if (fillRight) B[i * (N + 1) + j] = F(i * h1, j * h2) + (2 / h1) * psi_R(j * h2);

            } else if (j == 0) { // нижняя граница
                double a = k((i + 0.5) * h1, j * h2) * w_x(w, i, j) - k((i - 0.5) * h1, j * h2) * w_x_back(w, i, j);
                double b = k(i * h1, (1.0 - 0.5) * h2) * w_y_back(w, i, 1);
                A[i * (N + 1) + j] = -(2 / h2) * b + q(i * h1, j * h2) * w[i * (N + 1) + j] - a / h1;
                if (fillRight) B[i * (N + 1) + j] = F(i * h1, j * h2) + (2 / h2) * psi_B(i * h1);

            } else if (j == N) { // верхняя граница
                double a = k((i + 0.5) * h1, j * h2) * w_x(w, i, j) - k((i - 0.5) * h1, j * h2) * w_x_back(w, i, j);
                double b = k(i * h1, (j - 0.5) * h2) * w_y_back(w, i, j);
                A[i * (N + 1) + j] = (2 / h2) * b + q(i * h1, j * h2) * w[i * (N + 1) + j] - a / h1;
                if (fillRight) B[i * (N + 1) + j] = F(i * h1, j * h2) + (2 / h2) * psi_T(i * h1);

            } else { // внутренние точки
                double a = k((i + 0.5) * h1, j * h2) * w_x(w, i, j) - k((i - 0.5) * h1, j * h2) * w_x_back(w, i, j);
                double b = k(i * h1, (j + 0.5) * h2) * w_y(w, i, j) - k(i * h1, (j - 0.5) * h2) * w_y_back(w, i, j);
                A[i * (N + 1) + j] = - a / h1 - b / h2 + q(i * h1, j * h2) * w[i * (N + 1) + j];
                if (fillRight) B[i * (N + 1) + j] = F(i * h1, j * h2);
            }
        }
    }
}

double scal(double* a, double* b) {
    double res = 0.0;
    for (int i = 0; i <= M; ++i) {
        double rho1 = (i == 0 || i == M) ? 0.5 : 1.0;
        for (int j = 0; j <= N; ++j) {
            double rho2 = (j == 0 || j == N) ? 0.5 : 1.0;
            res += rho1 * rho2 * a[i * (N + 1) + j] * b[i * (N + 1) + j];
        }
    }

    return h1 * h2 * res;
}

void sub(double* res, double* a, double* b, int count) {
    for (int i = 0; i < count; ++i) {
        res[i] = a[i] - b[i];
    }
}

int main(int argc, char* argv[])
{
    M = atoi(argv[1]);
    N = atoi(argv[2]);
    double eps = atof(argv[3]);

    h1 = 4.0 / M;
    h2 = 3.0 / N;

    int dots_count = (M + 1) * (N + 1);
    double* Aw = new double[dots_count];
    double* Ar = new double[dots_count];
    double* w = new double[dots_count];
    double* next_w = new double[dots_count];
    double* r = new double[dots_count];
    double* B = new double[dots_count];
    double* diff = new double[dots_count];

    for (int i = 0; i <= M; ++i) {
        for (int j = 0; j <= N; ++j) {
            w[i * (N + 1) + j] = 0.0;
            next_w[i * (N + 1) + j] = 0.0;
        }
    }

    fill(Aw, w, B, true);

    double max1;
    do {
        for (int i = 0; i < dots_count; ++i) {
            w[i] = next_w[i];
        }
        fill(Aw, w, B);
        sub(r, Aw, B, dots_count);
        fill(Ar, r, B);
        double tmp1 = scal(Ar, r);
        double tmp2 = scal(Ar, Ar);
        double tau = tmp1 / tmp2;
        for (int i = 0; i < dots_count; ++i) {
            next_w[i] = w[i] - tau * r[i];
        }
        sub(diff, next_w, w, dots_count);
        max1 = 0.0;
        for (int i = 0; i < dots_count; ++i) {
            if (abs(diff[i]) > max1) max1 = abs(diff[i]);
        }

    } while (max1 > eps);

    //ofstream out("C:/Users/andre/Desktop/test0.csv");

    for (int i = 0; i <= M; ++i) {
        for (int j = 0; j <= N; ++j) {
            //out << next_w[i * (N + 1) + j];
            //if (j != N) out << ",";
            cout << next_w[i * (N + 1) + j];
            if (j != N) cout << ",";
        }
        //out << '\n';
        cout << '\n';
    }

    double* exact_solution = new double[dots_count];
    for (int i = 0; i <= M; ++i) {
        for (int j = 0; j <= N; ++j) {
            exact_solution[i * (N + 1) + j] = u(i * h1, j * h2);
        }
    }

    sub(diff, next_w, exact_solution, dots_count);
    cout << "norm = " << sqrt(scal(diff, diff)) << endl;

    return 0;
}

