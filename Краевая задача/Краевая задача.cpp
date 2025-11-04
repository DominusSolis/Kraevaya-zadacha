#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>

using namespace std;

const int N = 1000;
const double ax = 0.0, bx = 1.0;
// Аналитическое решение
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

// Первая производная аналитического решения
double y_an_derivative(double x) {
    return exp(x) - 0.2 * exp(-4.0 * x) + 0.8 * x * exp(-4.0 * x) - (1.0 / 6.0) + (1.0 / 36.0) * exp(-x);
}

// Вторая производная аналитического решения
double y_an_second_derivative(double x) {
    return exp(x) + 0.8 * exp(-4.0 * x) - 3.2 * x * exp(-4.0 * x) + 0.8 * exp(-4.0 * x) - (1.0 / 36.0) * exp(-x);
}

void saveToFile(const vector<double>& x, const vector<double>& y, const vector<double>& y_analytical,
    const vector<double>& first_order_approx, const vector<double>& second_order_approx,
    double avg_error, double max_error, double first_order_avg_error, double second_order_avg_error) {
    ofstream outFile("results.txt");

    if (!outFile.is_open()) {
        cerr << "Ошибка: не удалось создать файл results.txt" << endl;
        return;
    }

    outFile << fixed << setprecision(6);
    outFile << "Результаты решения краевой задачи:" << endl;
    outFile << "x\t\tЧисленное решение\tАналитическое решение\tАппроксимация 1-го порядка\tАппроксимация 2-го порядка" << endl;
    outFile << "--------------------------------------------------------------------------------------------------------" << endl;

    for (int i = 0; i <= N; i++) {
        outFile << setw(6) << setprecision(3) << x[i] << "\t\t"
            << setw(15) << setprecision(8) << y[i] << "\t\t"
            << setw(15) << setprecision(8) << y_analytical[i] << "\t\t"
            << setw(15) << setprecision(8) << first_order_approx[i] << "\t\t"
            << setw(15) << setprecision(8) << second_order_approx[i] << endl;
    }

    outFile << endl << "СТАТИСТИКА ПОГРЕШНОСТЕЙ:" << endl;
    outFile << "Средняя абсолютная погрешность: " << avg_error << endl;
    outFile << "Максимальная абсолютная погрешность: " << max_error << endl;
    outFile << "Средняя погрешность аппроксимации 1-го порядка: " << first_order_avg_error << endl;
    outFile << "Средняя погрешность аппроксимации 2-го порядка: " << second_order_avg_error << endl;

    outFile.close();

    cout << "Результаты сохранены в файл results.txt" << endl;
}

int main() {
    setlocale(LC_ALL, "Russian");
    vector<double> x(N + 1), y(N + 1), A(N + 1), B(N + 1), C(N + 1), F(N + 1), aa(N + 1), bb(N + 1);
    vector<double> y_analytical(N + 1);
    vector<double> first_order_approx(N + 1), second_order_approx(N + 1);
    vector<double> first_order_error(N + 1), second_order_error(N + 1);
    double h, error, max_error = 0.0;
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

    B[0] = 1.0;
    A[0] = 0.0;
    F[0] = 35.0 / 36.0;

    B[N] = 1.0 + 1.0 / h;
    C[N] = -1.0 / h;
    F[N] = 2.0 * exp(1) - (1.0 / 6.0) * exp(-1) + (2.0 / 5.0) * exp(-4.0);

    aa[0] = -A[0] / B[0];
    bb[0] = F[0] / B[0];

    for (i = 1; i <= N; i++) {
        aa[i] = -A[i] / (C[i] * aa[i - 1] + B[i]);
        bb[i] = (F[i] - C[i] * bb[i - 1]) / (C[i] * aa[i - 1] + B[i]);
    }

    y[N] = (F[N] - bb[N - 1] * C[N]) / (B[N] + aa[N - 1] * C[N]);

    for (i = N - 1; i >= 0; i--) {
        y[i] = aa[i] * y[i + 1] + bb[i];
    }

    // Вычисление погрешностей и аппроксимаций
    error = 0;
    double first_order_avg_error = 0.0, second_order_avg_error = 0.0;
    int first_order_count = 0, second_order_count = 0;

    for (i = 0; i <= N; i++) {
        double abs_error = abs(y[i] - y_analytical[i]);
        error += abs_error;

        if (abs_error > max_error) {
            max_error = abs_error;
        }

        // Аппроксимация первого порядка (разностная производная)
        if (i < N) {
            // Производная вперед
            first_order_approx[i] = (y[i + 1] - y[i]) / h;
            double analytical_derivative = y_an_derivative(x[i]);
            first_order_error[i] = abs(first_order_approx[i] - analytical_derivative);
            first_order_avg_error += first_order_error[i];
            first_order_count++;
        }
        else {
            // Для последней точки используем разностную производную назад
            first_order_approx[i] = (y[i] - y[i - 1]) / h;
            double analytical_derivative = y_an_derivative(x[i]);
            first_order_error[i] = abs(first_order_approx[i] - analytical_derivative);
            first_order_avg_error += first_order_error[i];
            first_order_count++;
        }

        // Аппроксимация второго порядка (вторая производная) - ИСПРАВЛЕННЫЙ БЛОК
        if (i == 0) {
            // Левая граница: разность вперед второго порядка
            second_order_approx[i] = (y[i + 2] - 2 * y[i + 1] + y[i]) / (h * h);
        }
        else if (i == N) {
            // Правая граница: разность назад второго порядка  
            second_order_approx[i] = (y[i] - 2 * y[i - 1] + y[i - 2]) / (h * h);
        }
        else {
            // Внутренние точки: центральная разность
            second_order_approx[i] = (y[i + 1] - 2 * y[i] + y[i - 1]) / (h * h);
        }

        double analytical_second_derivative = y_an_second_derivative(x[i]);
        second_order_error[i] = abs(second_order_approx[i] - analytical_second_derivative);
        second_order_avg_error += second_order_error[i];
        second_order_count++;
    }

    error = error / (N + 1);
    first_order_avg_error = first_order_avg_error / first_order_count;
    second_order_avg_error = second_order_avg_error / second_order_count;

    // Вывод
    cout << fixed << setprecision(6);
    cout << "РЕЗУЛЬТАТЫ РЕШЕНИЯ КРАЕВОЙ ЗАДАЧИ:" << endl;
    cout << "===================================" << endl;

    // Вывод первых 10 и последних 10 точек для наглядности
    cout << "Первые 10 точек:" << endl;
    cout << "x\t\tЧисленное\tАналитическое\tАппр.1-го пор.\tАппр.2-го пор." << endl;
    cout << "----------------------------------------------------------------" << endl;
    for (i = 0; i < 10; i++) {
        cout << setw(6) << setprecision(3) << x[i] << "\t\t"
            << setw(10) << setprecision(6) << y[i] << "\t"
            << setw(10) << setprecision(6) << y_analytical[i] << "\t"
            << setw(10) << setprecision(6) << first_order_approx[i] << "\t"
            << setw(10) << setprecision(6) << second_order_approx[i] << endl;
    }

    cout << endl << "Последние 10 точек:" << endl;
    cout << "x\t\tЧисленное\tАналитическое\tАппр.1-го пор.\tАппр.2-го пор." << endl;
    cout << "----------------------------------------------------------------" << endl;
    for (i = N - 9; i <= N; i++) {
        cout << setw(6) << setprecision(3) << x[i] << "\t\t"
            << setw(10) << setprecision(6) << y[i] << "\t"
            << setw(10) << setprecision(6) << y_analytical[i] << "\t"
            << setw(10) << setprecision(6) << first_order_approx[i] << "\t"
            << setw(10) << setprecision(6) << second_order_approx[i] << endl;
    }

    // Вывод статистики
    cout << endl << "СТАТИСТИКА ПОГРЕШНОСТЕЙ:" << endl;
    cout << "===================================" << endl;
    cout << "Средняя абсолютная погрешность: " << error << endl;
    cout << "Максимальная абсолютная погрешность: " << max_error << endl;
    cout << "Средняя погрешность аппроксимации 1-го порядка: " << first_order_avg_error << endl;
    cout << "Средняя погрешность аппроксимации 2-го порядка: " << second_order_avg_error << endl;
    cout << "Шаг сетки h: " << h << endl;

    // Проверка аппроксимации дифференциального уравнения
    cout << endl << "ПРОВЕРКА АППРОКСИМАЦИИ ДУ:" << endl;
    cout << "===================================" << endl;

    double max_residual = 0.0;
    for (i = 1; i < N; i++) {
        // Вычисляем невязку: y'' + p(x)y' + q(x)y - f(x)
        double y_second = (y[i + 1] - 2 * y[i] + y[i - 1]) / (h * h);
        double y_first = (y[i + 1] - y[i - 1]) / (2 * h);
        double residual = y_second + p(x[i]) * y_first + q(x[i]) * y[i] - f1(x[i]);

        if (abs(residual) > max_residual) {
            max_residual = abs(residual);
        }
    }
    cout << "Максимальная невязка ДУ: " << max_residual << endl;

    // Сохранение результатов в файл
    saveToFile(x, y, y_analytical, first_order_approx, second_order_approx,
        error, max_error, first_order_avg_error, second_order_avg_error);

    return 0;
}
