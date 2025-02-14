#include "opt_alg.h"
#include<iomanip>

void lab0();

void lab1();

void lab2();

void lab3();

void lab4();

void lab5();

void lab6();

int main() {
    try {
        lab4();
    } catch (string EX_INFO) {
        cerr << "ERROR:\n";
        cerr << EX_INFO << endl << endl;
    }
    system("pause");
    return 0;
}

void lab0() {
    //Funkcja testowa
    double epsilon = 1e-2;
    int Nmax = 10000;
    matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
    solution opt;
    a(0) = -1;
    a(1) = 2;
    opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
    cout << opt << endl << endl;
    solution::clear_calls();

    //Wahadlo
    Nmax = 1000;
    epsilon = 1e-2;
    lb = 0;
    ub = 5;
    double teta_opt = 1;
    opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
    cout << opt << endl << endl;
    solution::clear_calls();

    //Zapis symulacji do pliku csv
    matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{m2d(opt.x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
    ofstream Sout("symulacja_lab0.csv");
    Sout << hcat(Y[0], Y[1]);
    Sout.close();
    Y[0].~matrix();
    Y[1].~matrix();
}

void lab1() {
    try {
        // OPTYMALIZACJE

        string FILE_PATH = R"(C:\Users\Ala\Desktop\optyma\Optimalization.csv)";

        fstream csvFile;
        const int numOfElements = 100;
        double minVal = -15, maxVal = 15;
        set<double> uniqueNums;

        default_random_engine randomEngine(random_device{}());
        uniform_real_distribution<double> realDistribution(minVal, maxVal);

        while (uniqueNums.size() < numOfElements) {
            uniqueNums.insert(realDistribution(randomEngine));
        }

        for (double x0: uniqueNums) {
            // Expansion
            double d = 4.25, alpha = 1.00043;
            int Nmax = 1000;

            unique_ptr<double[]> p(expansion(ff1T, x0, d, alpha, Nmax));

            if (!p) {
                cout << "Expansion returned null for x0: " << x0 << endl;
                continue;
            }

            try {
                csvFile.open(FILE_PATH, ios::app);
                if (!csvFile.is_open()) {
                    throw runtime_error("Could not open file");
                }
                csvFile << x0 << "," << p[0] << "," << p[1] << "," << solution::f_calls << ",";
                csvFile.close();
            } catch (const exception &e) {
                cout << "Error writing expansion results: " << e.what() << endl;
                csvFile.close();
                continue;
            }
            solution::clear_calls();

            // Fibonacci
            try {
                solution Fibonacci = fib(ff1T, p[0], p[1], 0.0001);

                csvFile.open(FILE_PATH, ios::app);
                if (!csvFile.is_open()) {
                    throw runtime_error("Could not open file");
                }
                csvFile << m2d(Fibonacci.x) << "," << m2d(Fibonacci.y) << "," << solution::f_calls << ",";
                csvFile.close();
            } catch (const exception &e) {
                cout << "Error in Fibonacci: " << e.what() << endl;
                csvFile.close();
                continue;
            }
            solution::clear_calls();

            // Lagrange
            try {
                solution Lagrange = lag(ff1T, p[0], p[1], 0.0001, 1e-09, Nmax);

                csvFile.open(FILE_PATH, ios::app);
                if (!csvFile.is_open()) {
                    throw runtime_error("Could not open file");
                }
                csvFile << m2d(Lagrange.x) << "," << m2d(Lagrange.y) << "," << solution::f_calls << "\n";
                csvFile.close();
            } catch (const exception &e) {
                cout << "Error in Lagrange: " << e.what() << endl;
                csvFile.close();
                continue;
            }
            solution::clear_calls();
        }

        // SYMULACJA
        double DA0 = 0.001, d = 0.001, alpha = 2.257, epsilon = 1e-10, gamma = 1e-10;
        int Nmax = 10000;

        unique_ptr<double[]> range(expansion(ff1R, DA0, d, alpha, Nmax));

        if (!range) {
            throw runtime_error("Expansion for simulation returned null");
        }

        solution::clear_calls();

        solution FibonacciSolution = fib(ff1R, range[0], range[1], epsilon);
        cout << "Optimal DA value for Fibonacci's method: " << FibonacciSolution.x << "\n";
        cout << "Value of result function: " << FibonacciSolution.y << "\n";
        cout << "Number of calls of result function: " << solution::f_calls << "\n\n";
        solution::clear_calls();

        solution LagrangeSolution = lag(ff1R, range[0], range[1], epsilon, gamma, Nmax);
        cout << "Optimal DA value for Lagrange's method: " << LagrangeSolution.x << "\n";
        cout << "Value of result function: " << LagrangeSolution.y << "\n";
        cout << "Number of calls of result function: " << solution::f_calls << "\n";
        solution::clear_calls();

        double Pa = 0.5, Va0 = 5, Vb0 = 1, Tb0 = 20, t_0 = 0, t_step = 1, t_end = 2000;
        matrix Y0 = matrix(3, 1, Va0);
        Y0(1) = Vb0;
        Y0(2) = Tb0;

        double Da = m2d(FibonacciSolution.x);
        unique_ptr<matrix[]> FibonacciSimulation(solve_ode(df1R, t_0, t_step, t_end, Y0, Da, Pa));

        FILE_PATH =
                R"(C:\Users\Dell\Desktop\stuff\studia\semestr 5\Optymalizacja\project\data\data1\FibonacciSimulation.csv)";
        csvFile.open(FILE_PATH, ios::app);
        if (!csvFile.is_open()) {
            throw runtime_error("Could not open file for Fibonacci simulation");
        }
        csvFile << FibonacciSimulation[1] << "\n";
        csvFile.close();

        solution::clear_calls();

        Da = m2d(LagrangeSolution.x);
        unique_ptr<matrix[]> LagrangeSimulation(solve_ode(df1R, t_0, t_step, t_end, Y0, Da, Pa));

        FILE_PATH =
                R"(C:\Users\Dell\Desktop\stuff\studia\semestr 5\Optymalizacja\project\data\data1\LagrangeSimulation.csv)";
        csvFile.open(FILE_PATH, ios::app);
        if (!csvFile.is_open()) {
            throw runtime_error("Could not open file for Lagrange simulation");
        }
        csvFile << LagrangeSimulation[1] << "\n";
        csvFile.close();

        solution::clear_calls();
    } catch (const exception &e) {
        cout << "Fatal error in lab1: " << e.what() << endl;
    } catch (...) {
        cout << "Unknown error in lab1" << endl;
    }
}


void lab2() {
    try {
        std::string FILE_PATH = R"(C:\Users\Ala\Desktop\optyma\Optimalization.csv)";
        std::ofstream csvFile(FILE_PATH, std::ios::out);
        if (!csvFile.is_open()) throw std::runtime_error("Could not open file");

        double minVal = -1.0, maxVal = 1.0;
        std::default_random_engine randomEngine(std::random_device{}());
        std::uniform_real_distribution<double> realDistribution(minVal, maxVal);

        double epsilon = 1e-6;
        int Nmax = 1000;
        matrix ud1, ud2;
        std::vector<double> step_sizes = {0.44, 0.22, 0.11};


            matrix x0(2, 1);
            x0(0, 0) = realDistribution(randomEngine);
            x0(1, 0) = realDistribution(randomEngine);


                solution Xopt_HJ = HJ(ff2T, x0, 0.22, 0.45, epsilon, Nmax, ud1, ud2);
                solution::clear_calls();

                matrix s0(2, 1, 0.22);
                solution Xopt_Rosen = Rosen(ff2T, x0, s0, 1.22, 0.55, epsilon, Nmax, ud1, ud2);
                solution::clear_calls();



        csvFile.close();
    } catch (const std::exception &e) {
        std::cerr << "Fatal error in lab2: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown error in lab2" << std::endl;
    }
}

void lab3() {
    //FUNKCJA TESTOWA
         double penalty_start = 1.23;
         double penalty_scale_ext = 2.0;
         double penalty_scale_int = 1.5;
         double epsilon = 1e-3;
         int Nmax = 10000;
         int num_iterations = 100;
         std::vector<double> a_values = {4.0, 4.4934, 5.0};

         std::ofstream results("C:/Users/Ala/Desktop/optyma/results.csv");
         results << std::fixed << std::setprecision(7);
         results << "a;x_1;x_2;x1_ext;x2_ext;r_ext;y_ext;calls_ext;x1_int;x2_int;r_int;y_int;calls_int\n";

         for (double a_val : a_values) {
             matrix a_value = matrix(1, 1, a_val);

             for (int i = 0; i < num_iterations; ++i) {
                 matrix x0 = rand_mat(2, 1);
                 x0(0, 0) = x0(0, 0) * 2.0 - 1.0;
                 x0(1, 0) = x0(1, 0) * 2.0 - 1.0;

                 //Optymalizacja zewn�trzna
                 solution::clear_calls();
                 solution external_opt = pen(ff3T, x0, penalty_start, penalty_scale_ext, epsilon, Nmax, a_value,
                                              matrix(2, new double[2]{1.0, 2.0}));
                 double x1_ext = external_opt.x(0, 0);
                 double x2_ext = external_opt.x(1, 0);
                 double r_ext = norm(external_opt.x);
                 double y_ext = external_opt.y(0, 0);
                 int calls_ext = solution::f_calls;

                 //Optymalizacja wewn�trzna
                 solution::clear_calls();
                 solution internal_opt = pen(ff3T, x0, penalty_start, penalty_scale_int, epsilon, Nmax, a_value,
                                              matrix(2, new double[2]{1.0, 0.5}));
                 double x1_int = internal_opt.x(0, 0);
                 double x2_int = internal_opt.x(1, 0);
                 double r_int = norm(internal_opt.x);
                 double y_int = internal_opt.y(0, 0);
                 int calls_int = solution::f_calls;

                 // Zapis wynik�w w jednej linii
                 results << a_val << ";"
                         <<x0(0, 0) << ";" <<x0(1, 0) << ";"
                         << x1_ext << ";" << x2_ext << ";" << r_ext << ";" << y_ext << ";" << calls_ext << ";"
                         << x1_int << ";" << x2_int << ";" << r_int << ";" << y_int << ";" << calls_int << "\n";
             }
         }

    results.close();
    std::cout << "Optymalizacja zako�czona, wyniki zapisano w pliku results.csv" << std::endl;

    //SYMULACJA PROBLEMU
    matrix x_1 = matrix(2, 1);
    double c0 = 1.23;
    x_1(0) = 0.;     //vx_0  [-10, 10] (m/s)
    x_1(1) = 0.;     //omega [-15, 15] (rad/s)

    x_1 = pen(ff3R, x_1, c0, 2, epsilon, Nmax).x;
    cout<<endl;

    ff3R(x_1, c0, 2);
}

void lab4() {

    // Parametry optymalizacji
    matrix x0(2, 1, 0.0); // Punkt startowy [0, 0]^T
    x0(0, 0) = -10;       // Ustawienie x1 = -10
    x0(1, 0) = 10;        // Ustawienie x2 = 10

    double epsilon = 1e-6; // Dok�adno��
    int Nmax = 1000;       // Maksymalna liczba iteracji
    double h0 = 1.0;       // Pocz�tkowa d�ugo�� kroku

    matrix ud1, ud2;       // Dodatkowe dane u�ytkownika (tu puste)

    // Optymalizacja metod� Newtona
    cout << "Optymalizacja metod� Newtona:" << endl;
    solution result_newton = Newton(ff4T, gf4T, nullptr, x0, h0, epsilon, Nmax, ud1, ud2);
    cout << "Wynik Newton: x = " << result_newton.x
         << ", f(x) = " << result_newton.y << endl;

    // Optymalizacja metod� SD (najszybszego spadku)
    cout << "Optymalizacja metod� SD:" << endl;
    solution result_sd = SD(ff4T, gf4T, x0, h0, epsilon, Nmax, ud1, ud2);
    cout << "Wynik SD: x = " << result_sd.x
         << ", f(x) = " << result_sd.y << endl;

    // Optymalizacja metod� CG (gradient�w sprz�onych)
    cout << "Optymalizacja metod� CG:" << endl;
    solution result_cg = CG(ff4T, gf4T, x0, h0, epsilon, Nmax, ud1, ud2);
    cout << "Wynik CG: x = " << result_cg.x
         << ", f(x) = " << result_cg.y << endl;


}

void lab5() {
}

void lab6() {
}
