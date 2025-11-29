#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

template <typename Func>
double trapezium_rule(Func f, double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; ++i) {
        sum += f(a + i * h);
    }
    return sum * h;
}

double px, py, std_dev;
int num_sub_ints;

double integrate_triangle(double t1x, double t1y, double t2x, double t2y, double t3x, double t3y) {
    vector<vector<double>> M_inv = {
        {t2x - t1x, t3x - t2x},
        {t2y - t1y, t3y - t2y}
    };

    // Compute inverse of M_inv (2x2 matrix)
    double det = M_inv[0][0] * M_inv[1][1] - M_inv[0][1] * M_inv[1][0];
    if (abs(det) < 1e-12) {
        // cerr << "M_inv doesn't have an inverse\n";
        return 0;
    }

    vector<vector<double>> M = {
        { M_inv[1][1] / det, -M_inv[0][1] / det },
        { -M_inv[1][0] / det, M_inv[0][0] / det }
    };

    vector<double> t1 = {t1x, t1y};
    vector<double> p = {px, py};

    // b = -M * t1
    // commented out because it doesn't end up being relevant for calculation
    // vector<double> b = {
    //     - (M[0][0] * t1[0] + M[0][1] * t1[1]),
    //     - (M[1][0] * t1[0] + M[1][1] * t1[1])
    // };

    double CoefOfTransformation = abs(det);

    double A = M_inv[0][0];
    double B = M_inv[0][1];
    double C = t1[0] - p[0];
    double D = M_inv[1][0];
    double E = M_inv[1][1];
    double F = t1[1] - p[1];

    double T = B*B + E*E;

    double G, H;
    if (T == 0) {
        G = 1;
        H = 1;
    } else {
        G = (A*B + D*E) / T;
        H = (B*C + E*F) / T;
    }

    double L = C*C + F*F - T * H * H;
    double J = A*A + D*D - T * G * G;
    double K = 2 * (A*C + D*F - T*G*H);

    // define integrand as a lambda function
    auto integrand = [&](double u) {
        double exponent = -(J*u*u + K*u + L) / (2 * std_dev * std_dev);
        double var = sqrt(T) / (sqrt(2) * std_dev);
        double erf1_arg = var * (u + G*u + H);
        double erf2_arg = var * (G*u + H);
        return exp(exponent) * (erf(erf1_arg) - erf(erf2_arg));
    };

    double integral;
    if (num_sub_ints <= 0) {
        // Use a library integrator
        // REPLACE:
        integral=1;
    } else {
        // do trapezium rule
        integral = trapezium_rule(integrand, 0.0, 1.0, num_sub_ints);
    }

    return (CoefOfTransformation / (2 * sqrt(2 * M_PI * T) * std_dev)) * integral;
}

int main(int argc, char *argv[]) {
    // Usage:
    // [progname] <path_to_file_with_triangles> <px> <py> <std_dev> <num_sub_ints>
    px = atof(argv[2]);
    py = atof(argv[3]);
    std_dev = atof(argv[4]);
    num_sub_ints = atoi(argv[5]);

    // Read from the text file
    ifstream infile(argv[1]);

    if (!infile) {
        cerr << "Failed to open file.\n";
        return 1;
    }

    double sum = 0;

    string line;
    while (getline(infile, line)) {
        double t1x, t1y, t2x, t2y, t3x, t3y;

        // Parse the first pair
        infile >> t1x >> t1y;
        // Skip the '|'
        infile.ignore(1);
        // Parse the second pair
        infile >> t2x >> t2y;
        // Skip the '|'
        infile.ignore(1);
        // Parse the third pair
        infile >> t3x >> t3y;

        sum += integrate_triangle(t1x, t1y, t2x, t2y, t3x, t3y);
    }

    cout.precision(10);
    cout << sum << '\n';
}
