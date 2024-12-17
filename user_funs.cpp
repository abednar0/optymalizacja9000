#include"user_funs.h"
#include <functional>

matrix ff0T(matrix x, matrix ud1, matrix ud2) {
    matrix y(1, 1);
    y(0) = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
    return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{m2d(x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
    int n = get_len(Y[0]);
    double teta_max = Y[1](0, 0);
    for (int i = 1; i < n; ++i)
        if (teta_max < Y[1](i, 0))
            teta_max = Y[1](i, 0);
    y = abs(teta_max - m2d(ud1));
    Y[0].~matrix();
    Y[1].~matrix();
    return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(2, 1);
    double m = 1, l = 0.5, b = 0.5, g = 9.81;
    double I = m * pow(l, 2);
    dY(0) = Y(1);
    dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;
    return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
    return {-cos(0.1 * m2d(x)) * exp(-pow((0.1 * m2d(x) - 2 * M_PI), 2)) + 0.002 * pow(0.1 * m2d(x), 2)};
}

matrix ff1R(matrix x, matrix ud1, matrix ud2) {
    matrix Y0(3, 1);
    Y0(0) = 5.0;
    Y0(1) = 1.0;
    Y0(2) = 20.0;
    double t0 = 0.0;
    double dt = 1.0;
    int timesteps = 2000;

    matrix *Y = solve_ode(df1R, t0, dt, timesteps, Y0, x(0), x);

    double Tmax = Y[1](0, 2);
    for (int i = 1; i <= timesteps; ++i) {
        if (Y[1](i, 2) > Tmax) {
            Tmax = Y[1](i, 2);
        }
    }
    delete[] Y;

    return {fabs(Tmax - 50.0)};
}

matrix df1R(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix derivatives(3, 1);
    double a = 0.98, b = 0.63, g = 9.81, PA = 0.5, PB = 1, DB = 0.00365665;
    double F_in = 0.01, T_in = 20, TA = 90.0;
    double DA = m2d(ud1(0));

    double FA_out, FB_out;

    if (Y(0) > 0)
        FA_out = a * b * DA * sqrt(2 * g * Y(0) / PA);
    else
        FA_out = 0;

    if (Y(1) > 0)
        FB_out = a * b * DB * sqrt(2 * g * Y(1) / PB);
    else
        FB_out = 0;

    derivatives(0) = -FA_out;
    derivatives(1) = FA_out + F_in - FB_out;
    derivatives(2) = FA_out / Y(1) * (TA - Y(2)) + F_in / Y(1) * (T_in - Y(2));

    return derivatives;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2) {
    matrix y(1, 1);
    y(0) = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
    return y;
}

matrix df2R(double t, matrix Y, matrix ud1, matrix ud2) {
    double alpha = Y(0);
    double omega = Y(1);

    double k1 = ud2(0);
    double k2 = ud2(1);

    double m_arm = 1.0;
    double m_weight = 5.0;
    double l = 1.0;
    double g = 9.81;
    double b = 0.5;

    double I = (m_arm + m_weight) * pow(l, 2);

    double M = k1 * (M_PI - alpha) + k2 * (0.0 - omega);

    matrix dY(2, 1);
    dY(0) = omega;
    dY(1) = (M - m_arm * g * l * sin(alpha) - b * omega) / I;

    return dY;
}


matrix ff2R(matrix x, matrix ud1, matrix ud2) {
    double k1 = x(0);
    double k2 = x(1);

    double m = 1.0, l = 1.0, g = 9.81, b = 0.5;
    double alpha_target = M_PI;
    double omega_target = 0.0;

    matrix Y0(2, 1);
    Y0(0) = 0.0;
    Y0(1) = 0.0;

    matrix control_params(2, 1);
    control_params(0) = k1;
    control_params(1) = k2;

    matrix *result = solve_ode(df2R, 0.0, 0.1, 100.0, Y0, matrix(), control_params);

    double Q = 0.0;
    for (int i = 0; i < get_len(result[0]); ++i) {
        double alpha = result[1](i, 0);
        double omega = result[1](i, 1);
        double M = k1 * (alpha_target - alpha) + k2 * (omega_target - omega);
        Q += 10 * pow(alpha_target - alpha, 2) + pow(omega_target - omega, 2) + pow(M, 2);
    }

    delete[] result;

    return {Q};
}

// bool g1(matrix x1) {
//     return (-x1() + 1 <= 0);
// }
//
// bool g2(matrix x2) {
//     return (-x2() + 1 <= 0);
// }
//
// bool g3(matrix x1, matrix x2, double alpha) {
//     return (sqrt(pow(x1(), 2) + pow(x2(), 2)) - alpha <= 0);
// }

// matrix fR3(matrix x, matrix ud1, matrix ud2) {
//     matrix y;
//     matrix Y0(4, new double[4]{ 0, x(0), 100, 0 });
//     matrix* Y = solve_ode(fd3T, 0, 0.01, 7, Y0, ud1, x(1));
//     int n = get_len(Y[0]);
//     int i0 = 0, i50 = 0;
//
//     for (int i = 0; i < n; i++) {
//         if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50))
//             i50 = i;
//         if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2)))
//             i0 = i;
//     }
//
//     y = -Y[1](i0, 0);
//
//     if (abs(x(0)) - 10 > 0)
//         y = y + ud2 * pow(abs(x(0)) - 10, 2);
//     if (abs(x(1)) - 25 > 0)
//         y = y + ud2 * pow(abs(x(1)) - 25, 2);
//     if (abs(Y[1](koniec, 0) - 5) - 1 > 0)
//         y = y + ud2 * pow(abs(Y[1](koniec, 0) - 5) - 1, 2);
//
//     return y;
// }
//
// matrix ff3T(matrix x, matrix ud1, matrix ud2) {
//     double argument = M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2));
//     matrix y = sin(argument) / argument;
//
//     // Wektor ogranicze�
//     std::vector<std::function<double(matrix)>> constraints = {
//         [](matrix x) { return -x(0) + 1; },           // g1: -x1 + 1 <= 0
//         [](matrix x) { return -x(1) + 1; },           // g2: -x2 + 1 <= 0
//         [&ud1](matrix x) { return norm(x) - ud1(0); } // g3: ||x|| - a <= 0
//     };
//
//     // Zewn�trzna kara
//     if (ud2(1) > 1) {
//         for (const auto &g : constraints) {
//             double violation = std::max(0.0, g(x));
//             y = y + ud2(0) * pow(violation, 2);
//         }
//     }
//     // Wewn�trzna kara
//     else {
//         for (const auto &g : constraints) {
//             double violation = g(x);
//             if (violation >= 0) {
//                 return matrix(1, 1, 1e10); // Punkt poza obszarem dopuszczalnym
//             }
//             y = y - ud2(0) / (-violation);
//         }
//     }
//
//     return y;
// }
//
// matrix ff3R(matrix x, matrix ud1, matrix ud2) {
//
//     // Parametry wej�ciowe
//     double v_x = x(0, 0);    // Pr�dko�� pozioma [m/s]
//     double omega = x(1, 0);  // Pr�dko�� k�towa [rad/s]
//
//     // Sta�e fizyczne
//     const double m = 0.6;         // Masa pi�ki [kg]
//     const double r = 0.12;        // Promie� pi�ki [m]
//     const double rho = 1.2;       // G�sto�� powietrza [kg/m^3]
//     const double g = 9.81;        // Przyspieszenie grawitacyjne [m/s^2]
//     const double C = 0.47;        // Wsp�czynnik oporu
//     const double S = M_PI * r * r; // Powierzchnia przekroju poprzecznego
//
//     // Warunki pocz�tkowe
//     double y0 = 100.0;  // Pocz�tkowa wysoko�� [m]
//     double t0 = 0.0, dt = 0.01, t_max = 7.0;
//
//     // Pozycje i pr�dko�ci
//     double x_pos = 0.0, y_pos = y0;
//     double v_y = 0.0;  // Pocz�tkowa pr�dko�� pionowa
//     double t = t0;
//
//     // Symulacja ruchu
//     while (t <= t_max) {
//         // Pr�dko�� w powietrzu
//         double v_abs = sqrt(v_x * v_x + v_y * v_y);
//
//         // Si�y dzia�aj�ce na pi�k�
//         double D_x = 0.5 * C * rho * S * v_x * fabs(v_x);
//         double D_y = 0.5 * C * rho * S * v_y * fabs(v_y);
//         double F_x = rho * v_y * omega * M_PI * pow(r, 3);
//         double F_y = rho * v_x * omega * M_PI * pow(r, 3);
//
//         // Przyspieszenia
//         double a_x = -(D_x + F_x) / m;
//         double a_y = -(D_y + F_y + m * g) / m;
//
//         // Aktualizacja pr�dko�ci i pozycji
//         v_x += a_x * dt;
//         v_y += a_y * dt;
//         x_pos += v_x * dt;
//         y_pos += v_y * dt;
//
//         // Sprawdzenie czy pi�ka osi�gn�a ziemi�
//         if (y_pos <= 0) break;
//
//         t += dt;
//     }
//
//     // Funkcja celu: maksymalizacja x_pos
//     double goal = -x_pos;
//
//     // Penalizacja: Trafienie do kosza
//     double penalty = 0.0;
//
//     // Sprawdzanie odleg�o�ci od kosza (x = 5.0, y = 50.0)
//     double distance_to_hoop = sqrt(pow(x_pos - 5.0, 2) + pow(y_pos - 50.0, 2));
//     penalty += pow(distance_to_hoop > 0.5 ? distance_to_hoop - 0.5 : 0.0, 2);
//
//     // Dodatkowa kara za wyj�cie poza zakresy
//     if (v_x < -10 || v_x > 10) {
//         penalty += pow(v_x < -10 ? -10 - v_x : v_x - 10, 2);
//     }
//     if (omega < -15 || omega > 15) {
//         penalty += pow(omega < -15 ? -15 - omega : omega - 15, 2);
//     }
//
//     // Zwracamy wynik jako kombinacj� celu i kary
//     return matrix(goal + ud2(0, 0) * penalty);
//
// }

matrix df3R(double t, matrix Y, matrix ud1, matrix ud2) {
    double C = 0.47;
    double r = 0.12;
    double m = 0.6;
    double ro = 1.2;
    double g = 9.81;

    double S = M_PI * r * r;

    double omega = m2d(ud1);

    double x = Y(0);
    double y = Y(1);
    double vx = Y(2);
    double vy = Y(3);

    // Si�y
    double Dx = 0.5 * C * ro * S * vx * abs(vx);
    double Dy = 0.5 * C * ro * S * vy * abs(vy);
    double Fx = ro * vy * omega * M_PI * r * r * r;
    double Fy = ro * vx * omega * M_PI * r * r * r;

    matrix dY(4, 1);
    dY(0) = vx;
    dY(1) = vy;
    dY(2) = -(Dx + Fx) / m;
    dY(3) = -(Dy + Fy + m * g) / m;

    return dY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    matrix Y_temp(4, 1);
    Y_temp(0) = 0;
    Y_temp(1) = 100;
    Y_temp(2) = x(0);
    Y_temp(3) = 0;
    matrix* Y = solve_ode(df3R, 0, 0.01, 7, Y_temp, x(1));
    int n = get_len(Y[0]);
    int poczatek = 0, koniec = 0;

    for (int i = 0; i < n; i++) {
        if (abs(Y[1](i, 1) - 50) < abs(Y[1](koniec, 1) - 50))
            koniec = i;
        if (abs(Y[1](i, 1)) < abs(Y[1](poczatek, 1)))
            poczatek = i;
    }

    cout << Y[1](koniec, 0) << ";\t" << Y[1](koniec, 1) <<";\n";
    cout << Y[1](poczatek, 0) << ";\t" << Y[1](poczatek, 1) <<";\n";

    y = -Y[1](poczatek, 0);
    if (abs(x(0)) - 10 > 0)
        y = y + ud2 * pow(abs(x(0)) - 10, 2);
    if (abs(x(1)) - 15 > 0)
        y = y + ud2 * pow(abs(x(1)) - 15, 2);
    if (abs(Y[1](koniec, 0) - 5) > 0)
        y = y + ud2 * pow(abs(Y[1](koniec, 0) - 5), 2);

    return y;
}

matrix ff3T(matrix x, matrix ud1, matrix ud2) {
    double temp = M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2));
    matrix y = sin(temp) / temp;

    //Zewn�trzna
    if (ud2(1) > 1) {
        if (-x(0) + 1 > 0)
            y = y + (ud2)(0) * pow(-x(0) + 1, 2);
        if (-x(1) + 1 > 0)
            y = y + (ud2)(0) * pow(-x(1) + 1, 2);
        if (norm(x) - (ud1)(0) > 0)
            y = y + (ud2)(0) * pow(norm(x) - (ud1)(0), 2);
    }
    //Wewn�trzna
    else {
        if (-x(0) + 1 > 0)
            y = 1e12;
        else
            y = y - (ud2)(0) / (-x(0) + 1);

        if (-x(1) + 1 > 0)
            y = 1e12;
        else
            y = y - (ud2)(0) / (-x(1) + 1);

        if (norm(x) - (ud1)(0) > 0)
            y = 1e12;
        else
            y = y - (ud2)(0) / (sqrt(pow(x(0), 2) + pow(x(1), 2)) - (ud1)(0));
    }
    return y;
}

matrix ff4T(matrix x, matrix ud1, matrix ud2) {

        // Zak�adamy, �e x to wektor pionowy (2x1)
        double x1 = x(0, 0); // Pierwszy element wektora
        double x2 = x(1, 0); // Drugi element wektora

        // Funkcja celu
        double f = pow((x1 + 2 * x2 - 7), 2) + pow((2 * x1 + x2 - 5), 2);

        // Zwracamy wynik jako macierz 1x1 (zgodnie z konwencj� klasy `matrix`)
        return matrix(f);
}

matrix gf4T(matrix x, matrix ud1, matrix ud2) {
    try {
        double x1 = x(0, 0);
        double x2 = x(1, 0);

        double g1 = 2 * (x1 + 2 * x2 - 7) + 4 * (2 * x1 + x2 - 5);
        double g2 = 4 * (x1 + 2 * x2 - 7) + 2 * (2 * x1 + x2 - 5);

        matrix grad(2, 1);
        grad(0, 0) = g1;
        grad(1, 0) = g2;

        return grad;
    } catch (string ex_info) {
        throw ("matrix gf4T(...):\n" + ex_info);
    }
}

matrix fun4(matrix x, matrix ud1, matrix ud2) {
    return pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
}

matrix grad4(matrix x, matrix ud1, matrix ud2) {
    matrix g(2, 1);
    g(0) = 10 * x(0) + 8 * x(1) - 34;
    g(1) = 8 * x(0) + 10 * x(1) - 38;
    return g;
}

matrix hesj4(matrix x, matrix ud1, matrix ud2) {
    matrix H(2, 2);
    H(0, 0) = 10;
    H(0, 1) = 8;
    H(1, 0) = 8;
    H(1, 1) = 10;
    return H;
}

matrix fT4(matrix x, matrix ud1, matrix ud2) {
    matrix y;

    if (isnan(ud2(0, 0)))
        y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
    else
        y = fT4(ud2[0] + x * ud2[1], 0, ud1);
    return y;
}

// matrix fR4(matrix x, matrix ud1, matrix ud2) {
//     matrix y;
//     int m = 100;
//     int n = get_len(x);
//     static matrix X(n, m), Y(1, m);
//     if(solution::f_calls == 1) {
//         ifstream in("XData.txt");
//         in >> X;
//         in.close();
//         in.open("YData.txt");
//         in >> Y;
//         in.close();
//     }
//     int P = 0;
//     double h;
//     y = 0;
//     for (int i = 0; i < m; i++) {
//         h = m2d(trans(x) * X[i]);
//         h = 1.0 / (1.0 + exp(-h));
//         y = y - Y(0, i) * log(h) - (1 - Y(0, i)) * log(1 - h);
//     }
//     y = y / m;
//     return y;
// }
//
// matrix gf(matrix x, matrix ud1, matrix ud2) {
//     int m = 100;
//     int n = get_len(x);
//     matrix g(n, 1);
//     static matrix X(n, m), Y(1, m);
//     if (solution::g_calls == 1) {
//         ifstream in("XData.txt");
//         in >> X;
//         in.close();
//         in.open("YData.txt");
//         in >> Y;
//         in.close();
//     }
//
//     double h;
//     for (int j = 0; j < n; ++j) {
//         for (int i = 0; i < m; ++i) {
//             h = m2d(trans(x) * X[i]);
//             h = 1 / (1 + exp(-h));
//             g(j) = g(j) + X(j, i) * (h - Y(0, i));
//         }
//         g(j) = g(j) / m;
//     }
//     return g;
// }