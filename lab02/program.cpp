#include <iostream>
#include <map>
#include <string>
#include <chrono>
#include <iomanip>
#include <vector>

using namespace std;

const map<string, array<double, 3> > materials = {
    {"платина", {21450, 132.6, 69.7}},
    {"золото", {19320, 128.7, 317}},
    {"серебро", {10493, 235.4, 160}}
};



void simulate(string material, double L, double T0, double Tn, double time){
    double tau_values[] = {0.1, 0.01, 0.001, 0.0001};
    double h_values[] = {0.1, 0.01, 0.001, 0.0001};

    double rho   = materials.at(material)[0];
    double c     = materials.at(material)[1];
    double lambda = materials.at(material)[2];

    double temp_result[4][4] = {{0}};
    double time_result[4][4] = {{0}};

    for (int it = 0; it < 4; ++it) {
        double tau = tau_values[it];
        for (int ih = 0; ih < 4; ++ih) {
            double h = h_values[ih];
            auto start = chrono::steady_clock::now();

            int N = static_cast<int>(L / h); 

            vector<double> T_prev(N+1, T0);
            vector<double> A(N+1);
            vector<double> B(N+1);
            vector<double> C(N+1);
            vector<double> F(N+1);

            vector<double> alpha(N+1);
            vector<double> beta(N+1);

            double lambda_h2 = lambda / (h * h);
            double rho_c_tau = (rho * c) / tau;

            int nodes = static_cast<int>(L / h) + 1;
            int time_steps_total = static_cast<int>(time / tau); 

            for(size_t i = 0;i < N+1;i++) {
                A[i] = lambda_h2;
                C[i] = lambda_h2;
                B[i] = 2.0 * lambda_h2 + rho_c_tau;
            }

            for(size_t step = 0;step < time_steps_total;step++) {

                for(int i = 0;i < N+1;i++) {
                F[i] = -rho_c_tau*T_prev[i];
                }

                alpha[0] = 0.0;
                beta[0] = Tn;

                for(size_t i = 1;i < N;i++) {
                alpha[i] = A[i] / (B[i] - (C[i] * alpha[i-1]));
                beta[i] = (C[i]*beta[i-1] - F[i]) / (B[i] - (C[i] * alpha[i-1]));
                }

                vector<double> T_new(N+1);
                T_new[N] = Tn;

                for(int i = N - 1;i >= 0;i--) {
                T_new[i] = alpha[i] * T_new[i+1] + beta[i];
                }

                swap(T_prev,T_new);

            }

            auto end = chrono::steady_clock::now();

            chrono::duration<double> sim_time = end - start;

            temp_result[it][ih] = T_prev[N/2];
            time_result[it][ih] = sim_time.count();

            cout <<  "Расчёт: tau=" << fixed << setprecision(4) << tau << ", h=" << fixed << setprecision(4) << h << "...";
            cout << "T=" << fixed << setprecision(4) << temp_result[it][ih] << ", время=" << fixed << setprecision(4) << time_result[it][ih] << " с";
            cout << endl;
        }
    } 
    
    cout << "\nТемпература в центре пластины через " << time << " с:\n";
    cout << "tau\\h\t0.1\t0.01\t0.001\t0.0001\n";
    for (int i = 0; i < 4; ++i) {
        cout << tau_values[i] << "\t";
        for (int j = 0; j < 4; ++j) {
            cout << fixed << setprecision(4) << temp_result[i][j] << "\t";
        }
        cout << endl;
    }

    cout << "\nВремя расчёта (сек):\n";
    cout << "tau\\h\t0.1\t0.01\t0.001\t0.0001\n";
    for (int i = 0; i < 4; ++i) {
        cout << tau_values[i] << "\t";
        for (int j = 0; j < 4; ++j) {
            cout << fixed << setprecision(4) << time_result[i][j] << "\t";
        }
        cout << endl;
    }
}

int main(){
    double L, T0, Tn, time;
    string material;

    cout << "Введите материал пластины (платина, золото, серебро): ";
    cin >> material;
    cout << "Введите толщину пластины: ";
    cin >> L;
    cout << "Введите температуру внутри: ";
    cin >> T0;
    cout << "Введите температуру снаружи: ";
    cin >> Tn;
    cout << "Введите время моделирования: ";
    cin >> time;

    simulate(material, L, T0, Tn, time);
}