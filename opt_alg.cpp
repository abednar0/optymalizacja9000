#include "opt_alg.h"

solution
MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        while (true) {
            Xopt = rand_mat(N);
            for (int i = 0; i < N; ++i)
                Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
            Xopt.fit_fun(ff, ud1, ud2);
            if (Xopt.y < epsilon) {
                Xopt.flag = 1;
                break;
            }
            if (solution::f_calls > Nmax) {
                Xopt.flag = 0;
                break;
            }
        }
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution MC(...):\n" + ex_info);
    }
}

double *
expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2) {
    try {
        auto *p = new double[2]{0, 0};
        int i = 0;

        solution X0(x0), X1(x0 + d);
        X0.fit_fun(ff, ud1, ud2);
        X1.fit_fun(ff, ud1, ud2);

        if (X0.y == X1.y) {
            p[0] = m2d(X0.x);
            p[1] = m2d(X1.x);
            return p;
        }

        if (X1.y > X0.y) {
            d *= -1;
            X1.x = X0.x + d;
            X1.fit_fun(ff, ud1, ud2);
            if (X1.y >= X0.y) {
                p[0] = m2d(X1.x);
                p[1] = m2d(X0.x - d);
                return p;
            }
        }

        solution X2;
        do {
            X2.x = X0.x + pow(alpha, i) * d;
            X2.fit_fun(ff, ud1, ud2);

            if (solution::f_calls > Nmax) {
                p[0] = m2d(X0.x);
                p[1] = m2d(X2.x);
                return p;
            }

            if (X2.y >= X1.y)
                break;

            X0 = X1;
            X1 = X2;
            i++;
        } while (true);

        if (d > 0) {
            p[0] = m2d(X0.x);
            p[1] = m2d(X2.x);
        } else {
            p[0] = m2d(X2.x);
            p[1] = m2d(X0.x);
        }

        return p;
    }
    catch (string ex_info) {
        throw ("double* expansion(...):\n" + ex_info);
    }
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        unsigned int k = 1;

        vector<double> fibonacci;
        fibonacci.push_back(1);
        fibonacci.push_back(1);

        while (fibonacci[k] < (b - a) / epsilon) {
            fibonacci.push_back(fibonacci[k - 1] + fibonacci[k]);
            k++;
        }

        double c = b - fibonacci[k - 1] / fibonacci[k] * (b - a);
        solution C(c);
        C.fit_fun(ff, ud1, ud2);

        double d = a + b - c;
        solution D(d);
        D.fit_fun(ff, ud1, ud2);

        for (int i = 0; i < k - 3; i++) {
            if (C.y < D.y)
                b = d;
            else
                a = c;

            c = b - fibonacci[k - i - 2] / fibonacci[k - i - 1] * (b - a);
            C = solution(c);
            C.fit_fun(ff, ud1, ud2);

            d = a + b - c;
            D = solution(d);
            D.fit_fun(ff, ud1, ud2);
        }
        Xopt = solution(c);
        Xopt.fit_fun(ff, ud1, ud2);
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution fib(...):\n" + ex_info);
    }

}

solution
lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1,
    matrix ud2) {
    try {
        solution Xopt;

        solution A(a), B(b);
        double c = (a + b) / 2;
        solution C(c);

        A.fit_fun(ff, ud1, ud2);
        B.fit_fun(ff, ud1, ud2);
        C.fit_fun(ff, ud1, ud2);

        double d_prev = INFINITY;
        int iter = 0;

        while (true) {
            double l = m2d(A.y * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y * (pow(C.x(0), 2) - pow(A.x(0), 2)) +
                           C.y * (pow(A.x(0), 2) - pow(B.x(0), 2)));
            double m = m2d(A.y * (m2d(B.x) - m2d(C.x)) + B.y * (m2d(C.x) - m2d(A.x)) + C.y * (m2d(A.x) - m2d(B.x)));

            if (abs(m) <= epsilon) {
                Xopt = C;
                Xopt.flag = -1;
                return Xopt;
            }

            double d = 0.5 * (l / m);

            if (d < min(m2d(A.x), m2d(B.x)) || d > max(m2d(A.x), m2d(B.x))) {
                Xopt.flag = -2;
                return Xopt;
            }

            solution D(d);
            D.fit_fun(ff, ud1, ud2);

            if (solution::f_calls > Nmax) {
                Xopt = D;
                Xopt.flag = 0;
                return Xopt;
            }

            if (abs(d - d_prev) < gamma) {
                Xopt = D;
                Xopt.flag = 1;
                return Xopt;
            }

            if (d < m2d(C.x)) {
                if (m2d(D.y) < m2d(C.y)) {
                    B = C;
                    C = D;
                } else {
                    A = D;
                }
            } else {
                if (m2d(D.y) < m2d(C.y)) {
                    A = C;
                    C = D;
                } else {
                    B = D;
                }
            }

            double new_interval = abs(m2d(B.x) - m2d(A.x));
            Xopt.ud.add_row(new_interval);

            if (new_interval < epsilon) {
                Xopt = C;
                Xopt.flag = 1;
                return Xopt;
            }

            d_prev = d;
            iter++;
        }
    } catch (...) {
        throw "Error in lag";
    }
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        // Otwieranie pliku do zapisu iteracji
        std::ofstream logFile("HJ_Iterations.csv", std::ios::out);
        logFile << "Iteration;X0;X1;Func_Value\n";

        solution Xopt(x0);
        Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
        solution xB = Xopt;
        int f_calls = 0;

        int iteration = 0;
        while (s > epsilon) {
            solution x_trial = HJ_trial(ff, xB, s, ud1, ud2);

            if (x_trial.y < xB.y) {
                do {
                    xB = x_trial;
                    matrix new_point = xB.x * 2 - Xopt.x;
                    x_trial = solution(new_point);
                    x_trial.y = x_trial.fit_fun(ff, ud1, ud2);
                    f_calls++;

                    // Zapis bie��cej iteracji wewn�trznej do pliku
                    logFile << iteration << ";" << xB.x(0, 0) << ";" << xB.x(1, 0) << ";" << xB.y << "\n";

                    if (f_calls > Nmax) throw std::runtime_error("Exceeded max function calls");

                } while (x_trial.y < xB.y);

                Xopt = xB;
            } else {
                s *= alpha;
            }

            // Zapis bie��cej iteracji zewn�trznej do pliku, nawet je�li x_trial nie zosta� zaktualizowany
            logFile << iteration << ";" << xB.x(0, 0) << ";" << xB.x(1, 0) << ";" << xB.y << "\n";
            iteration++;

            if (f_calls > Nmax) throw std::runtime_error("Exceeded max function calls");
        }

        logFile.close();  // Zamkni�cie pliku po zako�czeniu iteracji
        return Xopt;
    } catch (string ex_info) {
        throw ("solution HJ(...):\n" + ex_info);
    }
}


solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2) {
    try {
        solution x_trial = XB;
        int* size = get_size(XB.x);
        int n = size[0];
        delete[] size;

        for (int j = 0; j < n; ++j) {
            matrix step_forward = XB.x;
            step_forward(j, 0) += s;
            solution trial_fwd(step_forward);
            trial_fwd.y = trial_fwd.fit_fun(ff, ud1, ud2);

            if (trial_fwd.y < x_trial.y) {
                x_trial = trial_fwd;
            } else {
                matrix step_backward = XB.x;
                step_backward(j, 0) -= s;
                solution trial_bwd(step_backward);
                trial_bwd.y = trial_bwd.fit_fun(ff, ud1, ud2);

                if (trial_bwd.y < x_trial.y) {
                    x_trial = trial_bwd;
                }
            }
        }

        return x_trial;
    }
    catch (string ex_info) {
        throw ("solution HJ_trial(...):\n" + ex_info);
    }
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        // Otwieranie pliku do zapisu iteracji
        std::ofstream logFile("Rosen_Iterations.csv", std::ios::out);
        logFile << "Iteration;X0;X1;Func_Value\n";

        solution Xopt(x0);
        Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
        matrix x = x0;
        matrix s = s0;
        int iter_count = 0;
        int n = get_len(x0);

        matrix dir = ident_mat(n);

        while (iter_count < Nmax) {
            solution x_trial = Xopt;

            for (int j = 0; j < n; ++j) {
                matrix step = s(j) * dir[j];
                solution forward_trial(x + step);
                forward_trial.y = forward_trial.fit_fun(ff, ud1, ud2);

                solution backward_trial(x - step);
                backward_trial.y = backward_trial.fit_fun(ff, ud1, ud2);

                if (forward_trial.y < Xopt.y) {
                    Xopt = forward_trial;
                    x = Xopt.x;
                    s(j) *= alpha;
                } else if (backward_trial.y < Xopt.y) {
                    Xopt = backward_trial;
                    x = Xopt.x;
                    s(j) *= alpha;
                } else {
                    s(j) *= beta;
                }
            }

            // Zapis bie��cej iteracji do pliku
            logFile << iter_count << ";" << x(0, 0) << ";" << x(1, 0) << ";" << Xopt.y << "\n";
            iter_count++;

            if (norm(s) < epsilon) break;
        }

        logFile.close();  // Zamkni�cie pliku po zako�czeniu iteracji
        return Xopt;
    } catch (string ex_info) {
        throw ("solution Rosen(...):\n" + ex_info);
    }
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        //solution Xopt;
        double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
        solution X(x0), X1;
        matrix c0(2, new double[2]{ c,dc });
        while (true) {
            X1 = sym_NM(ff, X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c0);
            if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax) {
                X1.flag = 0;
                return X1;
            }
            X = X1;
            c0(0) = c0(0) * dc;
        }
    } catch (std::string ex_info) {
        throw ("solution pen(...):\n" + ex_info);
    }
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        int dim = get_len(x0);
        matrix identity = ident_mat(dim);
        int vertices_count = dim + 1;
        solution* simplex = new solution[vertices_count];

        simplex[0].x = x0;
        simplex[0].fit_fun(ff, ud1, ud2);

        for (int i = 1; i < vertices_count; ++i) {
            simplex[i].x = simplex[0].x + s * identity[i - 1];
            simplex[i].fit_fun(ff, ud1, ud2);
        }

        solution reflected, expanded, contracted;
        matrix centroid;
        int idx_min, idx_max;

        while (true) {
            idx_min = idx_max = 0;
            for (int i = 1; i < vertices_count; ++i) {
                if (simplex[i].y(0) < simplex[idx_min].y(0))
                    idx_min = i;
                if (simplex[i].y(0) > simplex[idx_max].y(0))
                    idx_max = i;
            }

            centroid = matrix(dim, 1);
            for (int i = 0; i < vertices_count; ++i) {
                if (i != idx_max) centroid = centroid + simplex[i].x;
            }
            centroid = centroid / dim;

            reflected.x = centroid + alpha * (centroid - simplex[idx_max].x);
            reflected.fit_fun(ff, ud1, ud2);

            if (reflected.y(0) < simplex[idx_max].y(0) && simplex[idx_min].y(0) <= reflected.y(0)) {
                simplex[idx_max] = reflected;
            } else if (reflected.y(0) < simplex[idx_min].y(0)) {
                expanded.x = centroid + gamma * (reflected.x - centroid);
                expanded.fit_fun(ff, ud1, ud2);

                simplex[idx_max] = (expanded.y(0) < reflected.y(0)) ? expanded : reflected;
            } else {
                contracted.x = centroid + beta * (simplex[idx_max].x - centroid);
                contracted.fit_fun(ff, ud1, ud2);

                if (contracted.y(0) < simplex[idx_max].y(0)) {
                    simplex[idx_max] = contracted;
                } else {
                    for (int i = 0; i < vertices_count; ++i) {
                        if (i != idx_min) {
                            simplex[i].x = delta * (simplex[i].x + simplex[idx_min].x);
                            simplex[i].fit_fun(ff, ud1, ud2);
                        }
                    }
                }
            }

            double max_distance = norm(simplex[idx_min].x - simplex[0].x);
            for (int i = 1; i < vertices_count; ++i) {
                double dist = norm(simplex[idx_min].x - simplex[i].x);
                if (dist > max_distance) max_distance = dist;
            }

            if (max_distance < epsilon) return simplex[idx_min];
        }
    } catch (std::string ex_info) {
        throw ("solution sym_NM(...):\n" + ex_info);
    }
}

solution
SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        Xopt.x = x0; // Punkt pocz�tkowy
        int iteration = 0;

        // G��wna p�tla iteracyjna
        while (iteration < Nmax) {
            // Obliczanie gradientu w bie��cym punkcie
            Xopt.g = Xopt.grad(gf);

            // Wyznaczanie kierunku d(i) = -grad(f)
            matrix d = -Xopt.g;

            // Sprawdzenie normy gradientu (warunek zako�czenia)
            if (norm(d) < epsilon) {
                Xopt.flag = 1; // Minimum osi�gni�te
                break;
            }

            // Wyznaczanie d�ugo�ci kroku h(i) metod� ekspansji i z�otego podzia�u
            double *exp_bounds = expansion(ff, m2d(Xopt.x(0, 0)), 1.0, 1.2, Nmax, ud1, ud2);
            double a = exp_bounds[0];
            double b = exp_bounds[1];

            solution golden_result = golden(ff, a, b, epsilon, Nmax, ud1, ud2);
            double h = m2d(golden_result.x(0, 0)); // Optymalna d�ugo�� kroku

            // Aktualizacja punktu x(i+1) = x(i) + h(i) * d(i)
            Xopt.x = Xopt.x + h * d;

            // Obliczanie warto�ci funkcji celu w nowym punkcie
            Xopt.y = Xopt.fit_fun(ff);

            // Warunek stopu ||x(i+1) - x(i)|| < epsilon
            if (norm(d * h) < epsilon) {
                Xopt.flag = 1; // Minimum osi�gni�te
                break;
            }

            iteration++;
        }

        // Je�li liczba iteracji przekroczy�a Nmax
        if (iteration >= Nmax) {
            Xopt.flag = -1; // Niepowodzenie w osi�gni�ciu minimum
        }

        return Xopt;
    } catch (string ex_info) {
        throw ("solution SD(...):\n" + ex_info);
    }
}

solution
CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        Xopt.x = x0; // Punkt pocz�tkowy
        int iteration = 0;

        // Zmienne pomocnicze do przechowywania poprzedniego gradientu i kierunku
        matrix g_prev, d_prev;

        // G��wna p�tla iteracyjna
        while (iteration < Nmax) {
            // Obliczanie gradientu w bie��cym punkcie
            Xopt.g = Xopt.grad(gf);

            // Wyznaczanie kierunku d(i)
            matrix d;
            if (iteration == 0) {
                // Pierwszy krok: d(0) = -grad(f)
                d = -Xopt.g;
            } else {
                // Obliczanie wsp�czynnika Beta
                double beta = (trans(Xopt.g) * Xopt.g)(0, 0) / (trans(g_prev) * g_prev)(0, 0);
                d = -Xopt.g + beta * d_prev;
            }

            // Sprawdzenie normy gradientu (warunek zako�czenia)
            if (norm(d) < epsilon) {
                Xopt.flag = 1; // Minimum osi�gni�te
                break;
            }

            // Wyznaczanie d�ugo�ci kroku h(i) metod� ekspansji i z�otego podzia�u
            double *exp_bounds = expansion(ff, m2d(Xopt.x(0, 0)), 1.0, 1.2, Nmax, ud1, ud2);
            double a = exp_bounds[0];
            double b = exp_bounds[1];

            solution golden_result = golden(ff, a, b, epsilon, Nmax, ud1, ud2);
            double h = m2d(golden_result.x(0, 0)); // Optymalna d�ugo�� kroku

            // Aktualizacja punktu x(i+1) = x(i) + h(i) * d(i)
            Xopt.x = Xopt.x + h * d;

            // Obliczanie warto�ci funkcji celu w nowym punkcie
            Xopt.y = Xopt.fit_fun(ff);

            // Przechowywanie bie��cych warto�ci gradientu i kierunku
            g_prev = Xopt.g;
            d_prev = d;

            // Warunek stopu ||x(i+1) - x(i)|| < epsilon
            if (norm(d * h) < epsilon) {
                Xopt.flag = 1; // Minimum osi�gni�te
                break;
            }

            iteration++;
        }

        // Je�li liczba iteracji przekroczy�a Nmax
        if (iteration >= Nmax) {
            Xopt.flag = -1; // Niepowodzenie w osi�gni�ciu minimum
        }

        return Xopt;
    } catch (string ex_info) {
        throw ("solution CG(...):\n" + ex_info);
    }
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
                matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1,
                matrix ud2) {
    try {
        solution Xopt; // Obiekt, kt�ry przechowuje wynik optymalizacji
        Xopt.x = x0;   // Ustawienie punktu pocz�tkowego
        int iteration = 0; // Licznik iteracji

        // G��wna p�tla iteracyjna
        while (iteration < Nmax) {
            // Obliczanie gradientu i hesjanu w bie��cym punkcie
            Xopt.g = Xopt.grad(gf);
            Xopt.H = Xopt.hess(Hf);

            // Wyznaczanie kierunku d(i) = -H(x)^(-1) * grad(f)
            matrix d = -inv(Xopt.H) * Xopt.g;

            // Sprawdzenie normy kierunku (warunek stopu dla ||d|| < epsilon)
            if (norm(d) < epsilon) {
                Xopt.flag = 1; // Minimum osi�gni�te
                break;
            }

            // Wyznaczanie d�ugo�ci kroku h(i) metod� ekspansji i z�otego podzia�u
            double *exp_bounds = expansion(ff, m2d(Xopt.x(0, 0)), 1.0, 1.2, Nmax, ud1, ud2);
            double a = exp_bounds[0];
            double b = exp_bounds[1];

            solution golden_result = golden(ff, a, b, epsilon, Nmax, ud1, ud2);
            double h = m2d(golden_result.x(0, 0)); // Optymalna d�ugo�� kroku

            // Aktualizacja punktu x(i+1) = x(i) + h(i) * d(i)
            Xopt.x = Xopt.x + h * d;

            // Obliczanie warto�ci funkcji celu w nowym punkcie
            Xopt.y = Xopt.fit_fun(ff);

            // Sprawdzenie warunku stopu ||x(i+1) - x(i)|| < epsilon
            if (norm(d * h) < epsilon) {
                Xopt.flag = 1; // Minimum osi�gni�te
                break;
            }

            iteration++;
        }

        // Je�li liczba iteracji przekroczy�a Nmax
        if (iteration >= Nmax) {
            Xopt.flag = -1; // Niepowodzenie w osi�gni�ciu minimum
        }

        return Xopt; // Zwracamy wynik optymalizacji
    } catch (string ex_info) {
        throw ("solution Newton(...):\n" + ex_info);
    }
}

solution
golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt; // Obiekt przechowuj�cy wynik optymalizacji
        const double alpha = (sqrt(5.0) - 1.0) / 2.0; // Sta�a z�otego podzia�u
        int iteration = 0;

        // Pocz�tkowe warto�ci przedzia�u i punkt�w
        double c = b - alpha * (b - a);
        double d = a + alpha * (b - a);

        // Obliczenie warto�ci funkcji celu w punktach c i d
        matrix xc = c;
        matrix xd = d;

        double fc = m2d(solution().fit_fun(ff));
        double fd = m2d(solution().fit_fun(ff));

        // G��wna p�tla iteracyjna
        while ((b - a > epsilon) && (iteration < Nmax)) {
            if (fc < fd) {
                // Przesuwamy przedzia� na lewo
                b = d;
                d = c;
                c = b - alpha * (b - a);

                // Przeliczenie warto�ci funkcji celu w nowym punkcie
                xd = xc;
                fd = fc;

                xc = c;
                fc = m2d(solution().fit_fun(ff));
            } else {
                // Przesuwamy przedzia� na prawo
                a = c;
                c = d;
                d = a + alpha * (b - a);

                // Przeliczenie warto�ci funkcji celu w nowym punkcie
                xc = xd;
                fc = fd;

                xd = d;
                fd = m2d(solution().fit_fun(ff));
            }

            iteration++;
        }

        // Wyznaczenie punktu minimalnego (�rodek przedzia�u [a, b])
        double x_opt = (a + b) / 2.0;
        Xopt.x = x_opt; // Wsp�rz�dna optymalna
        Xopt.y = solution().fit_fun(ff); // Warto�� funkcji w minimum
        Xopt.flag = (iteration < Nmax) ? 1 : -1; // Flaga sukcesu (1) lub niepowodzenia (-1)

        return Xopt;
    } catch (string ex_info) {
        throw ("solution golden(...):\n" + ex_info);
    }
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution Powell(...):\n" + ex_info);
    }
}

solution
EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution EA(...):\n" + ex_info);
    }
}
