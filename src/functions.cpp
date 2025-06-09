#include "functions.hpp"

extern "C" {
    void dgtsv_(const int* n, const int* nrhs, double* dl, double* d, double* du, double* b, const int* ldb, int* info);
}


void wczytaj_dane(double &D, double &d, double &a, double &T, int &ciecz, double &mi, double &delta_t, double &tk) {
    std::cout << "Podaj srednice rury: ";
    std::cin >> D;

    std::cout << "Podaj srednice rdzenia: ";
    std::cin >> d;

    std::cout << "Podaj amplitude sily: ";
    std::cin >> a;

    std::cout << "Podaj okres ruchu: ";
    std::cin >> T;

    std::cout << "Podaj rodzaj cieczy:\n"
              << "\t 1 - woda\n"
              << "\t 2 - rtec\n"
              << "\t 3 - powietrze\n"
              << "\t 4 - wodor\n";
    std::cin >> ciecz;

    switch (ciecz) {
        case 1:
            mi = 0.89 * 0.001;
            delta_t = pow(10, -5);
            break;
        case 2:
            mi = 1.554 * 0.001;
            delta_t = pow(10, -4);
            break;
        case 3:
            mi = 17.08 * 0.000001;
            delta_t = pow(10, -2);
            break;
        case 4:
            mi = 21.37 * 0.000001;
            delta_t = pow(10, -2);
            break;
        default:
            std::cerr << "Blad: nieznana ciecz!" << std::endl;
            exit(-1);
    }

    std::cout << "Podaj czas symulacji: ";
    std::cin >> tk;
}

void macierz_ukladu(double delta_t, double mi, double h, const std::vector<double>& r, int N,  std::vector<std::vector<double>> &A) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j)
			A[i][j] = 0;
	}
	A[0][0] = 1. / delta_t + 2. * mi / (h * h);
	A[0][1] = -mi * (r[0] + h / 2.) / (r[0] * h * h);

	for (int i = 1; i < N - 1; ++i) {
		A[i][i - 1] = -mi * (r[i] - h / 2.) / (r[i] * h * h);
		A[i][i] = 1. / delta_t + 2. * mi / (h * h);
		A[i][i + 1] = -mi * (r[i] + h / 2.) / (r[i] * h * h);
	}
	A[N - 1][N - 1 - 1] = -mi * (r[N - 1] - h / 2.) / (r[N - 1] * h * h);
	A[N - 1][N - 1] = 1. / delta_t + 2. * mi / (h * h);

}

void wektor_prawych_stron(double a, double T, double delta_t, double t, const std::vector<double>& u0, int N,  double *b) {
	double w;
	w = 2. * 4. * atan(1.0) / T;

	for (int i = 0; i < N; ++i) {
		b[i] = -a * sin(w * t) + 1. / delta_t * u0[i];
	}

}

void rozwiazanie_calkowite(const std::vector<double>& u1, int N,  std::vector<double>& u) {

	for (int i = 1; i < N - 1; ++i)
		u[i] = u1[i - 1];

}

void srednia_predkosc(const std::vector<double>& u, const std::vector<double>& r, int n,  double* u_sr, double h) {
	double u_calka_t;
	trapz(u, h, n, &u_calka_t);
	*u_sr = u_calka_t / (r[n - 1] - r[0]);
}

void przekroj_rury(double r, double R, double* A) {
	*A = 4 * atan(1.) * (R * R - r * r);
}

void wydatek(double u_sr, double A, double* Q) {
	*Q = A * u_sr;
}

void trapz(const std::vector<double>& x, double h, int n, double *y) {
    *y = 0;
    for (int i = 0; i < n-1; ++i) {
        *y += (x[i] + x[i + 1]) * h / 2;
    }
}

void Thomas(const std::vector<std::vector<double>> &A,  std::vector<double>& x,  double* b, int N) {
    double *alpha, *beta;
    alpha = (double*)malloc(N * sizeof(double));
    beta = (double*)malloc(N * sizeof(double));
    alpha[0] = -A[0][1] / A[0][0];
    beta[0] = b[0] /  A[0][0];
    for (int i = 0; i < N-2; ++i) {
        alpha[i+1] = -A[i + 1][i +1 + 1] / (A[i+1][i+1] + alpha[i]*A[i+1][i+1-1]);
        beta[i+1] = (b[i+1] -A[i+1][i+1-1] * beta[i]) / (A[i + 1][i + 1] + alpha[i] * A[i + 1][i + 1 - 1]);
    }
    x[N - 1] = (b[N - 1] - A[N - 1][N - 1 -1] * beta[ N - 2]) / (A[N - 1][N - 1]  + alpha[N-2]* A[N - 1][N - 1 - 1]);
    for (int i = N - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
    free(alpha);
    free(beta);
}

void extractTridiagonals(const std::vector<std::vector<double>> &A, 
                           std::vector<double> &dl, 
                           std::vector<double> &diag, 
                           std::vector<double> &du) {
   size_t n = A.size();
    if(n == 0) {
        throw std::invalid_argument("Macierz jest pusta!");
    }
    
    // Sprawdzenie, czy macierz jest kwadratowa.
    for (size_t i = 0; i < n; ++i) {
        if (A[i].size() != n) {
            throw std::invalid_argument("Macierz nie jest kwadratowa!");
        }
    }
    
    // Ustalenie rozmiarów wektorów wynikowych:
    // - diag: n elementów
    // - dl i du: n-1 elementów (dla n > 1)
    diag.resize(n);
    if (n > 1) {
        dl.resize(n - 1);
        du.resize(n - 1);
    } else {
        dl.clear();
        du.clear();
    }
    
    // Wypełnienie przekątnej:
    for (size_t i = 0; i < n; ++i) {
        diag[i] = A[i][i];
    }
    
    // Wypełnienie dolnej i górnej przekątnej:
    for (size_t i = 0; i < n - 1; ++i) {
        dl[i] = A[i + 1][i];   // elementy poniżej głównej przekątnej
        du[i] = A[i][i + 1];   // elementy powyżej głównej przekątnej
    }
}


void wizualizacja(sf::RenderWindow &window, const std::vector<double>& u, const std::vector<double>& r, double r1, double r2, int n, int ciecz) {
    double k = 1000.0;
    double x1 = 100.0;
    double y1 = 100.0;
    double x2 = 600.0;
    double y2 = 100.0;

    sf::Color lineColor = sf::Color::White;
    sf::Vertex line1[] = {
        sf::Vertex(sf::Vector2f(x1, y1), lineColor),
        sf::Vertex(sf::Vector2f(x2, y2), lineColor)
    };
    window.draw(line1, 2, sf::Lines);

    sf::Vertex line2[] = {
        sf::Vertex(sf::Vector2f(x1, y1 + 2 * r2 * k), lineColor),
        sf::Vertex(sf::Vector2f(x2, y2 + 2 * r2 * k), lineColor)
    };
    window.draw(line2, 2, sf::Lines);

    sf::Vertex line3[] = {
        sf::Vertex(sf::Vector2f(x1, y1 + r2 * k - r[0] * k), lineColor),
        sf::Vertex(sf::Vector2f(x2, y1 + r2 * k - r[0] * k), lineColor)
    };
    window.draw(line3, 2, sf::Lines);

    sf::Vertex line4[] = {
        sf::Vertex(sf::Vector2f(x1, y1 + r2 * k + r[0] * k), lineColor),
        sf::Vertex(sf::Vector2f(x2, y1 + r2 * k + r[0] * k), lineColor)
    };
    window.draw(line4, 2, sf::Lines);
    
   
    sf::Color circleColor(100, 0, 0);

    for (int i = 0; i < n; ++i) {
		
        float centerX = static_cast<float>((x1 + x2) / 2.0 + u[i]);
        float centerYUpper = static_cast<float>(y1 + r2 * k - r[i] * k);
        float centerYLower = static_cast<float>(y1 + r2 * k + r[i] * k);
        
      
        sf::CircleShape circle(2);

        circle.setOrigin(2, 2);
        circle.setFillColor(circleColor);

    
        circle.setPosition(centerX, centerYUpper);
        window.draw(circle);

   
        circle.setPosition(centerX, centerYLower);
        window.draw(circle);
    }
}

void wizualizacja_thread_func(sf::RenderWindow &window,
                              const std::vector<double>& u,
                              const std::vector<double>& r,
                               double r1,  double r2,  int n,  int ciecz){
    window.setActive(true);
    while (window.isOpen()&& simulationRunning.load()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        std::vector<double> u_copy;
        std::vector<double> r_copy;
        {
            std::lock_guard<std::mutex> lock(dataMutex);
            u_copy = u;
            r_copy = r;
        
		}
            window.clear(sf::Color::Black);
            wizualizacja(window, u_copy, r_copy, r1, r2, n, ciecz);
            window.display();   
    }
}

void symuluj_dgtsv(double a, double T,double delta_t, double t, std::vector<double>& u01, std::vector<double>& u1,std::vector<double>& u, 
                     const std::vector<double>& r_j,     int n1, int n,  double *b,const std::vector<double>&dl,
                     const std::vector<double>&diag,  const std::vector<double>&du,double u_sr, double h, double A_r, double Q,
                    FILE* fb,FILE* fu, FILE* fr,FILE* fu_srQ, int nt ){
    std::vector<double> dl_copy(n1 - 1);
    std::vector<double> diag_copy(n1);
    std::vector<double> du_copy(n1 - 1);

      const int aa = 1;
    int info = 0;
    for (unsigned int j = 0; j < nt; ++j) {
        sf::sleep(sf::milliseconds(5));
        std::lock_guard<std::mutex> lock(dataMutex);
		//for (int i = 0; i < static_cast<int>(u.size()); ++i)
           // std::cout << "b[" << i << "] = " << b[i] << std::endl;
        wektor_prawych_stron(a, T, delta_t, t, u01, n1, b);
		std::ranges::copy(dl,std::begin(dl_copy));
		std::ranges::copy(diag,std::begin(diag_copy));
		std::ranges::copy(du,std::begin(du_copy));
        dgtsv_(&n1, &aa, dl_copy.data(), diag_copy.data(), du_copy.data(), b, &n1, &info);
        for (size_t i = 0; i < static_cast<size_t>(n1); ++i) {
            u1.at(i) = b[i];
        }
		//Thomas(A, u1, b, n1);
        std::cout << "Czas = " << t << " s" << std::endl;
       // for (int i = 0; i < static_cast<int>(u.size()); ++i)
           // std::cout << "u[" << i << "] = " << u[i] << std::endl;

        rozwiazanie_calkowite(u1, n, u);
        srednia_predkosc(u, r_j, n, &u_sr, h);
        wydatek(u_sr, A_r, &Q);

        // Zapis do plikow
        if (j % 100 == 0 && t > 0.649) {
            for (int i = 0; i < n1; ++i) {
                fprintf(fb, "%lf\n", b[i]);
            }
            for (int i = 0; i < n; ++i) {
                fprintf(fu, "%lf\n", u[i]);
            }
            for (int i = 0; i < n; ++i) {
                fprintf(fr, "%lf\t%lf\n", r_j[i], t);
            }
            fprintf(fu_srQ, "%lf\t%lf\t%lf\n", t, u_sr, Q);
            printf("%d\n", j);
        }

        // Nadpisanie stanu z poprzedniej iteracji
        for (int i = 0; i < n1; ++i) {
            u01[i] = u1[i];
        }
        t += delta_t;
    }
}

void symuluj_Thomas(double a, double T,double delta_t, double t, std::vector<double>& u01, std::vector<double>& u1,std::vector<double>& u, 
                     const std::vector<double>& r_j, const std::vector<std::vector<double>>& A,   int n1, int n,  double *b,double u_sr, double h, double A_r, double Q,
                    FILE* fb,FILE* fu, FILE* fr,FILE* fu_srQ, int nt){

    for (unsigned int j = 0; j < nt; ++j) {
        sf::sleep(sf::milliseconds(5));
        std::lock_guard<std::mutex> lock(dataMutex);
		//for (int i = 0; i < static_cast<int>(u.size()); ++i)
           // std::cout << "b[" << i << "] = " << b[i] << std::endl;
        wektor_prawych_stron(a, T, delta_t, t, u01, n1, b);

		Thomas(A, u1, b, n1);
        std::cout << "Czas = " << t << " s" << std::endl;
       // for (int i = 0; i < static_cast<int>(u.size()); ++i)
           // std::cout << "u[" << i << "] = " << u[i] << std::endl;
        rozwiazanie_calkowite(u1, n, u);
        srednia_predkosc(u, r_j, n, &u_sr, h);
        wydatek(u_sr, A_r, &Q);
        // Zapis do plikow
        if (j % 100 == 0 && t > 0.649) {
            for (int i = 0; i < n1; ++i) {
                fprintf(fb, "%lf\n", b[i]);
            }
            for (int i = 0; i < n; ++i) {
                fprintf(fu, "%lf\n", u[i]);
            }
            for (int i = 0; i < n; ++i) {
                fprintf(fr, "%lf\t%lf\n", r_j[i], t);
            }
            fprintf(fu_srQ, "%lf\t%lf\t%lf\n", t, u_sr, Q);
            printf("%d\n", j);
        }

        // Nadpisanie stanu z poprzedniej iteracji
        for (int i = 0; i < n1; ++i) {
            u01[i] = u1[i];
        }
        t += delta_t;
    }
}