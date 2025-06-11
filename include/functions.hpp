#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <exception>
#include <iostream>
#include <thread>
#include <iomanip>
#include <SFML/Graphics.hpp>
#include "globals.hpp"
#include <algorithm>
#include <ranges>
#include <cblas.h>
  
      


void wczytaj_dane(double &D, double &d, double &a, double &T, int &ciecz, double &mi, double &delta_t, double &tk);

void macierz_ukladu(double delta_t, double mi, double h, const std::vector<double>& r, int N, std::vector<std::vector<double>>&A);

void wektor_prawych_stron(double a, double T, double delta_t, double t, const std::vector<double>& u0, int N, double* b);

void rozwiazanie_calkowite(const std::vector<double>& u, int N,  std::vector<double>& u_c);

void srednia_predkosc(const std::vector<double>& u,const std::vector<double>& r, int n,  double *u_sr, double h);

void przekroj_rury(double r, double R, double* A);

void wydatek(double u_sr, double A, double* Q);

void trapz(const std::vector<double>& x, double h, int n, double *y);

void Thomas(const std::vector<std::vector<double>> &A,  std::vector<double>& x,  double*  b, int N);

void Thomas_3D(const std::vector<double> &dl,const std::vector<double> &diag,const std::vector<double> &du,  std::vector<double>& x,  double* b, int N);


void extractTridiagonals(const std::vector<std::vector<double>> &A, std::vector<double> &dl,  std::vector<double> &diag, std::vector<double> &du);




void wizualizacja(sf::RenderWindow &window, const std::vector<double>&u, const std::vector<double>& r, double r1, double r2, int n, int ciecz);

void wizualizacja_thread_func(sf::RenderWindow &window,
                              const std::vector<double>& u,
                              const std::vector<double>& r,
                             double r1,  double r2,  int n,  int ciecz);


void symuluj_dgtsv(double a, double T,double delta_t, double t, std::vector<double>& u01, std::vector<double>& u1,std::vector<double>& u, 
                     const std::vector<double>& r_j,     int n1, int n,  double *b,const std::vector<double>&dl,
                     const std::vector<double>&diag,  const std::vector<double>&du,double u_sr, double h, double A_r, double Q,
                     int nt );

void symuluj_Thomas(double a, double T,double delta_t, double t, std::vector<double>& u01, std::vector<double>& u1,std::vector<double>& u, 
                     const std::vector<double>& r_j, const std::vector<std::vector<double>>& A, int n1, int n,  double *b,double u_sr, double h, double A_r, double Q,
                    int nt);

void symuluj_Thomas3d(double a, double T,double delta_t, double t, std::vector<double>& u01, std::vector<double>& u1,std::vector<double>& u, 
                     const std::vector<double>& r_j, const std::vector<double>&dl,
                     const std::vector<double>&diag,  const std::vector<double>&du,   int n1, int n,  double *b,double u_sr, double h, double A_r, double Q, int nt);

#endif 