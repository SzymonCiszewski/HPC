#include "functions.hpp"
#include "globals.hpp"
#include "benchmark/benchmark.h"


void symulacja(benchmark::State& state) {
    FILE* fA = fopen("macierz_ukladu.txt", "w");
    FILE* fb = fopen("wektor_prawych_stron.txt", "w");
    FILE* fu = fopen("wektor_niewiadomych.txt", "w");
    FILE* fr = fopen("wektor_r_wezlow.txt", "w");
    FILE* fu_srQ = fopen("czas_srednia_predkosc_wydatek.txt", "w");

    // Srednica rury [m]
    double D;
    // Srednica rdzenia [m]
    double d;
    // Krok pomiedzy wezlami [m]
    double h = 1e-3;
    // Pole przekroju rury [m^2]
    double A_r;
    // Amplituda sily [N]
    double a;
    // Rodzaj cieczy
    int ciecz;
    // Okres ruchu [s]
    double T;
    // Czas [s]
    double t = 0;
    // Krok czasowy [s]
    double delta_t;
    // Czas trwania symulacji [s]
    double tk;
    // Lepkosc dynamiczna cieczy [Pa*s]
    double mi;

    // Wczytanie danych od uzytkownika
    wczytaj_dane(D, d, a, T, ciecz, mi, delta_t, tk);

    // Promien rdzenia wewnetrznego [m]
    double r = d / 2;
    // Promien rury [m]
    double R = D / 2;
    przekroj_rury(r, R, &A_r);

    // Liczba krokow czasowych [-]
    int nt = (tk - t) / delta_t;
    // Liczba wezlow na y(r) – liczba wierszy macierzy wezlow [-]
    int n = (R - r) / h + 1;
    // Liczba wezlow z niewiadoma u [-] (bez warunkow brzegowych)
    const int n1 = n - 2;

    printf("Liczba wezlow, krokow czasowych : %d\t%d\n", n, nt);

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

    macierz_ukladu(delta_t, mi, h, r_j1, n1, A);

    std::vector<double> dl(n1 - 1);
    std::vector<double> diag(n1);
    std::vector<double> du(n1 - 1);


    try {
        extractTridiagonals(A, dl, diag, du);
        for (size_t i = 0; i < dl.size(); ++i)
            std::cout << "dl[" << i << "] = " << dl[i] << std::endl;
        for (size_t i = 0; i < diag.size(); ++i)
            std::cout << "diag[" << i << "] = " << diag[i] << std::endl;
        for (size_t i = 0; i < du.size(); ++i)
            std::cout << "du[" << i << "] = " << du[i] << std::endl;
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Wyjatek std::out_of_range: " << e.what() << std::endl;
    }
    catch (...) {
        std::cerr << "Rzucono nieznany wyjatek" << std::endl;
    }

  
	
    sf::RenderWindow window(sf::VideoMode(1000, 500), "Wizualizacja SFML");
    window.setFramerateLimit(120);
    window.setActive(false);
	
    std::thread vizThread(
        wizualizacja_thread_func,
        std::ref(window),
        std::cref(u),
        std::cref(r_j),
        r, R, n, ciecz
    );
	for(auto _ : state) // 3. Benchmark loop
		{
         // 4. Code to be benchmarked
        
           //symuluj_Thomas( a,  T, delta_t,  t,  u01,u1, u, r_j,  A,    n1,  n,   b, u_sr,  h,  A_r,  Q, fb, fu,  fr, fu_srQ,  nt);
		   symuluj_dgtsv(a, T, delta_t,  t, u01,  u1,u,r_j, n1,  n,   b,dl,diag, du, u_sr,  h,  A_r,  Q,fb, fu, fr, fu_srQ,  nt );
           // benchmark::DoNotOptimize(output);
            benchmark::ClobberMemory();
        }
         
     
    printf("Symulacja zakonczona\n");
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n1; ++j) {
            fprintf(fA, "%lf\t", A[i][j]);
        }
        fprintf(fA, "\n");
    }
	
    simulationRunning.store(false);
    if (vizThread.joinable())
        vizThread.join();

    free(b);
    fclose(fA);
    fclose(fb);
    fclose(fu);
    fclose(fr);
    fclose(fu_srQ);
    sf::sleep(sf::seconds(5));


}


BENCHMARK(symulacja);
BENCHMARK_MAIN();