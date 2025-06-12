#include "functions.hpp"
#include "globals.hpp"
#include "benchmark/benchmark.h"



void symulacja(benchmark::State& state) {
    
    // Liczba wezlow na y(r) – liczba wierszy macierzy wezlow [-]
    const int n=state.range(0);
    const int work_per_iteration =41;

    // Srednica rury [m]
    double D=0.2;
    // Srednica rdzenia [m]
    double d=0.02;
    // Pole przekroju rury [m^2]
    double A_r;
    // Amplituda sily [N]
    double a=50;
    // Rodzaj cieczy
    int ciecz;
    // Okres ruchu [s]
    double T=10;
    // Czas [s]
    double t = 0;
    // Krok czasowy [s]
    double delta_t=pow(10, -2);
    // Czas trwania symulacji [s]
    double tk=10;
    // Lepkosc dynamiczna cieczy [Pa*s]
    double mi=17.08 * 0.000001;
    // Promien rdzenia wewnetrznego [m]
    double r = d / 2;
    // Promien rury [m]
    double R = D / 2;
    przekroj_rury(r, R, &A_r);
    // Liczba krokow czasowych [-]
    int nt = (tk - t) / delta_t;
    // Liczba wezlow na y(r) – liczba wierszy macierzy wezlow [-]
    //int n;// = (R - r) / h + 1;
    // Liczba wezlow z niewiadoma u [-] (bez warunkow brzegowych)
    const int n1 = n - 2;
     // Krok pomiedzy wezlami [m]
    double h = (R-r)/(n-1);

    std::vector<double> u(n);     // Wektor u dla wszystkich wezlow
    std::vector<double> u1(n1);     // Wektor u dla wezlow bez scianek
    std::vector<double> u0(n);      // Wektor warunku poczatkowego
    std::vector<double> u01(n1);    // Wektor warunku poczatkowego bez scianek
    std::vector<double> r_j(n);     // Wektor wspolrzednych r dla wszystkich wezlow
    std::vector<double> r_j1(n1);   // Wektor wspolrzednych r dla wezlow bez scianek
    std::vector<std::vector<double>> A(n1, std::vector<double>(n1));

    double u_sr, Q;
    double* b = (double*)malloc(n1 * sizeof(double));

    for (int i = 0; i < n; ++i) {
        u0[i] = 0;
        r_j[i] = r + i * h;
    }
    for (int i = 0; i < n1; ++i) {
        u01[i] = 0;
        r_j1[i] = r + (i + 1) * h;
    }
    // Warunki brzegowe
    u[0] = 0;
    u[n - 1] = 0;
  
   // macierz_ukladu(delta_t, mi, h, r_j1, n1, A);
    std::vector<double> dl(n1 - 1);
   std::vector<double> diag(n1);
   std::vector<double> du(n1 - 1);
    macierz_ukladu_3d(delta_t,mi,h,r_j1,n1,dl,diag,du);

  
    /*
    sf::RenderWindow window(sf::VideoMode(1000, 500), "Wizualizacja SFML");
    window.setFramerateLimit(120);
    window.setActive(false);
	
    std::thread vizThread(
        wizualizacja_thread_func,
        std::ref(window),
        std::cref(u),
        std::cref(r_j),
        r, R, n, ciecz
    );*/
	for(auto _ : state) // 3. Benchmark loop
		{
         // 4. Code to be benchmarked
        
           //symuluj_Thomas( a,  T, delta_t,  t,  u01,u1, u, r_j,  A,    n1,  n,   b, u_sr,  h,  A_r,  Q,nt);
           //symuluj_Thomas3d( a,  T, delta_t,  t,  u01,u1, u, r_j,  dl,diag,du,    n1,  n,   b, u_sr,  h,  A_r,  Q,nt);
		   symuluj_dgtsv(a, T, delta_t,  t, u01,  u1,u,r_j, n1,  n,   b,dl,diag, du, u_sr,  h,  A_r,  Q, nt );
           // benchmark::DoNotOptimize(output);
            benchmark::ClobberMemory();
        }
         /*
    simulationRunning.store(false);
    if (vizThread.joinable())
        vizThread.join();*/

    free(b);
   
    //sf::sleep(sf::seconds(5));
    state.SetBytesProcessed(state.iterations() * work_per_iteration);
}


BENCHMARK(symulacja)->RangeMultiplier(2)->Range(1<<10, 1<<15);
BENCHMARK_MAIN();