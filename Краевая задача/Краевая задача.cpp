// Я УСТАЛ ЕГО ПИСАТЬ, ОН У МЕНЯ НЕ РАБОТАЛ, Я НЕ ПОНИМАЮ ПОЧЕМУ ТАКАЯ БОЛЬШАЯ ПОГРЕШНОСТЬ, МОЖЕТ БЫТЬ НЕПРАВИЛЬНО ВЫЧИСЛИЛ КОЭФФИЦИЕНТЫ...Т_Т
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>

using namespace std;

const int N = 1000;
const double ax = 0.0, bx = 1.0;

double y_an(double x) {
    return exp(x) - 0.2 * x * exp(-4.0 * x) - ((x / 6.0) + (1.0 / 36.0) * exp(-x));
}

double p(double x) {
    return 3.0;
}

double q(double x) {
    return -4.0;
}

double f1(double x) {
    return exp(-4.0 * x) + x * exp(-x);
}

void saveToFile(const vector<double>& x, const vector<double>& y, const vector<double>& y_analytical, double error) {
    ofstream outFile("results.txt");

    if (!outFile.is_open()) {
        cerr << "Ошибка: не удалось создать файл results.txt" << endl;
        return;
    }

    outFile << fixed << setprecision(6);
    outFile << "Результаты решения краевой задачи:" << endl;
    outFile << "x\t\tЧисленное решение\tАналитическое решение" << endl;
    outFile << "------------------------------------------------" << endl;

    for (int i = 0; i <= N; i++) {
        outFile << setw(6) << setprecision(2) << x[i] << "\t\t"
            << setw(10) << setprecision(6) << y[i] << "\t\t"
            << setw(10) << setprecision(6) << y_analytical[i] << endl;
    }

    outFile << endl << "Средняя ошибка: " << error << endl;
    outFile.close();

    cout << "Результаты сохранены в файл results.txt" << endl;
}

int main() {
    setlocale(LC_ALL, "Russian");
    vector<double> x(N + 1), y(N + 1), A(N + 1), B(N + 1), C(N + 1), F(N + 1), aa(N + 1), bb(N + 1);
    vector<double> y_analytical(N + 1);
    double h, error;
    int i;

    h = (bx - ax) / N;

    for (i = 0; i <= N; i++) {
        x[i] = ax + h * i;
        y_analytical[i] = y_an(x[i]);
    }

    for (i = 0; i <= N - 1; i++) {
        C[i] = 1.0 / (h * h) - p(x[i]) / (2 * h);
        A[i] = 1.0 / (h * h) + p(x[i]) / (2 * h);
        B[i] = -2.0 / (h * h) + q(x[i]);
        F[i] = f1(x[i]);
    }

    // Граничные условия
    B[0] = 1.0;
    A[0] = 0.0;
    F[0] = 35.0 / 36.0;

    B[N] = 1.0 + 1.0 / h;
    C[N] = -1.0 / h;
    F[N] = 2.0 * exp(1) - (1.0 / 6.0) * exp(-1) + (2.0 / 5.0) * exp(-4.0);
    // Прямой ход прогонки
    aa[0] = -A[0] / B[0];
    bb[0] = F[0] / B[0];

    for (i = 1; i <= N; i++) {
        aa[i] = -A[i] / (C[i] * aa[i - 1] + B[i]);
        bb[i] = (F[i] - C[i] * bb[i - 1]) / (C[i] * aa[i - 1] + B[i]);
    }

    // Обратный ход прогонки
    y[N] = (F[N] - bb[N - 1] * C[N]) / (B[N] + aa[N - 1] * C[N]);

    for (i = N - 1; i >= 0; i--) {
        y[i] = aa[i] * y[i + 1] + bb[i];
    }

    // Вывод результатов
    cout << fixed << setprecision(6);
    for (i = 0; i <= N; i++) {
        cout << setw(6) << setprecision(2) << x[i] << " , "
            << setw(10) << setprecision(6) << y[i] << " , "
            << setw(10) << setprecision(6) << y_analytical[i] << endl;
    }

    // Вычисление средней ошибки
    error = 0;
    for (i = 1; i <= N; i++) {
        error += abs(y[i] - y_analytical[i]);
    }
    error = error / N;
    cout << "Средняя ошибка: " << error << endl;

    // Сохранение результатов в файл
    saveToFile(x, y, y_analytical, error);

    return 0;
}