#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
using namespace std;

double f_R(double m, double a, double b, double f_a, double f_b) {
	return (m * (b - a) + ((f_b - f_a) * (f_b - f_a)) / (m * (b - a)) - 2 * (f_b + f_a));
}

struct Result {
	double x;
	double y;	
	int k;
	double m;
	double time = 0.0;
	Result(double x_, double y_, int k_, double m_) {
		x = x_;
		y = y_;		
		k = k_;
		m = m_;
	}
	const Result& operator=(const Result& obj) {
		if (this == &obj) return (*this);
		x = obj.x;		y = obj.y;
		k = obj.k;		m = obj.m;
		return (*this);
	}	friend ostream& operator<<(ostream& os, const Result& obj)
	{
		os << "arg_extremum = " << obj.x << "\nobj_func_extremum = " << obj.y
			<< "\ncount operation to search = " << obj.k << "\ntime = " << obj.time
			<< "\nConst L = " << obj.m << endl << endl;		return os;
	}
};

Result algorithm_for_searching(double a, double b, double (*f)(double x), double r, int max_count_operation, double eps) {
	int count_operation = 1;  // k - total points 
	vector<double>  arg = { a, b }, objective_function = { f(a), f(b) };  

	// creating the first two segments
	double max_M = abs((f(b) - f(a)) / (b - a)), m = 0;
	if (max_M == 0)
		m = 1;
	else
		m = r * max_M;
	double x_k = 0.5 * (b + a) - (objective_function[1] - objective_function[0]) / (2 * m);
	arg.insert(arg.cbegin() + 1, x_k);
	objective_function.insert(objective_function.cbegin() + 1, f(x_k));

	double	res_arg = arg[0], res_objective_function = objective_function[0];

	while (count_operation < max_count_operation) {
		int ind_max_R = 1;  // number x of the end of the segment with max R
		max_M = abs((objective_function[1] - objective_function[0]) / (arg[1] - arg[0]));
		for (int i = 2; i < arg.size(); i++) {
			double M = abs((objective_function[i] - objective_function[i - 1]) / (arg[i] - arg[i - 1]));
			if (M > max_M)
				max_M = M;

		}
		if (max_M == 0)
			m = 1;
		else
			m = r * max_M;
		// calculate max R
		double max_R = f_R(m, arg[0], arg[1], objective_function[0], objective_function[1]);
		for (int i = 2; i < arg.size(); i++) {
			double R = f_R(m, arg[i-1], arg[i], objective_function[i-1], objective_function[i]);
			if (max_R < R) {
				max_R = R;
				ind_max_R = i;
			}
		}

		x_k = 0.5 * (arg[ind_max_R] + arg[ind_max_R - 1]) - (objective_function[ind_max_R] - objective_function[ind_max_R - 1]) / (2 * m);
		arg.insert(arg.cbegin() + ind_max_R, x_k);
		objective_function.insert(objective_function.cbegin() + ind_max_R , f(x_k));
		if (res_objective_function > f(x_k)) {
			res_objective_function = f(x_k);
			res_arg = x_k;
		}
		count_operation++;

		if (arg[ind_max_R+1] - arg[ind_max_R - 1] < eps) {
			break;
		}
	}
	Result res = { res_arg, res_objective_function, count_operation, m };
	return res;
}

double f_sin(double x) {
	return -sin(x);
}

double f_x(double x) {
	return (x - 10) * (x - 10) + 3;
}

double f_1(double x) {
	return (sin(x) + sin(10 * x / 3));
}

double f_2(double x) {
	double res = 0;
	for (int k = 1; k < 6; k++) {
		res += k * sin((k + 1) * x + k);
	}
	return -res;
}

double f_3(double x) {
	return ((3 * x - 1.4) * sin(18 * x));
}

double f_4(double x) {
	return -((x + sin(x)) * exp(-(x * x)));
}

double f_5(double x) {
	return sin(x) + sin(10 * x / 3) + log(x) - 0.84 * x + 3;
}

double f_6(double x) {
	//const double PI = acos(-1.0);
	const double PI = 3.141592653589793;
	return -sin(2 * PI * x) * exp(-x);
}

double f_7(double x) {
	return (x * x - 5 * x + 6) / (x * x + 1);
}

double f_8(double x) {
	return (-x + sin(3 * x) - 1);
}

double f_9(double x) {
	return (2 * (x - 3) * (x - 3) + exp(x * x / 2));
}

int main() {
	double r = 1.5, error_of_arg = 0.001, max_count_operation = 100000;
	auto start = chrono::high_resolution_clock::now();
	Result res_f1 = algorithm_for_searching(2.7, 7.5, f_1, r, max_count_operation, error_of_arg);
	auto end = chrono::high_resolution_clock::now();
	double time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f1.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f2 = algorithm_for_searching(0.0, 10.0, f_2, r, max_count_operation, error_of_arg);
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f2.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f3 = algorithm_for_searching(0.0, 1.2, f_3, r, max_count_operation, error_of_arg);
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f3.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f4 = algorithm_for_searching(-10.0, 10.0, f_4, r, max_count_operation, error_of_arg);
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f4.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f5 = algorithm_for_searching(2.7, 7.5, f_5, r, max_count_operation, error_of_arg);
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f5.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f6 = algorithm_for_searching(0.0, 4.0, f_6, r, max_count_operation, error_of_arg); // not true
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f6.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f7 = algorithm_for_searching(-5.0, 5.0, f_7, r, max_count_operation, error_of_arg); // not true
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f7.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f8 = algorithm_for_searching(0.0, 6.5, f_8, r, max_count_operation, error_of_arg);
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f8.time = time_taken;

	start = chrono::high_resolution_clock::now();
	Result res_f9 = algorithm_for_searching(-3.0, 3.0, f_9, r, max_count_operation, error_of_arg);
	end = chrono::high_resolution_clock::now();
	time_taken =
		chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	time_taken *= 1e-9;
	res_f9.time = time_taken;

	cout << "f1:" << endl << res_f1
		<< "f2:" << endl << res_f2
		<< "f3:" << endl << res_f3
		<< "f4:" << endl << res_f4
		<< "f5:" << endl
		<< res_f5
		<< "f6:" << endl << res_f6
		<< "f7:" << endl << res_f7
		<< "f8:" << endl << res_f8
		<< "f9:" << endl << res_f9;
	return 0;
}
