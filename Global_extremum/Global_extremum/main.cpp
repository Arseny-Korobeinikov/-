
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;



double f_R(double m, double a, double b, double f_a, double f_b) {
	return (m * (b - a) + ((f_b - f_a) * (f_b - f_a)) / (m * (b - a)) - 2 * (f_b - f_a));
}


struct Result {
	double x;
	double y;
	int k;
	Result(double x_, double y_, int k_) {
		x = x_;
		y = y_;
		k = k_;
	}

	const Result& operator=(const Result& obj) {
		if (this == &obj) return (*this);
		x = obj.x;
		y = obj.y;
		k = obj.k;

		return (*this);
	}



	friend ostream& operator<<(ostream& os, const Result& obj)
	{
		os << "arg_extremum = " << obj.x << "\nobj_func_extremum = " << obj.y << "\ncount operation to search = " << obj.k << endl << endl;
		return os;
	}
};






Result algorithm_for_searching(double a, double b, double (*f)(double x), double r, int max_count_operation, double eps) {
	int count_operation = 1;  // k - total points 
	vector<double>  arg = { a, b }, objective_function = { f(a), f(b) };  //лучше назвать objective function

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



	Result res = { res_arg, res_objective_function, count_operation };
	return res;
}



double f_sin(double x) {
	return -sin(x);
}

double f_x(double x) {
	return (x -10)*(x-10)+3;
}


int main() {
	Result res = algorithm_for_searching(0, 3.1415, f_sin, 1.2, 100, 0.2);
	Result res1 = algorithm_for_searching(-2, 2, f_x, 1.2, 100, 0.2);

	cout << res << res1;

	return 0;
}		







/* как из этой версии сделать следующую версию
Сохранить и сделать следующую.
Придумать способ как ускорить вычисления
*/









		/*double M = abs((objective_function[t+i] - objective_function[t+i - 1]) / (arg[t+i] - arg[t+i - 1]));
			if (M > max_M) {
				max_M = M;
				flag_change_m = 1;
				if (max_M == 0)
					m = 1;
				else
					m = r * max_M;
			}*/





			//R.push_back( m * (x_k - a) + ( (objective_function[1] - objective_function[0]) * (objective_function[1] - objective_function[0])) / (m * (x_k - a)) - 2 * (objective_function[1] + objective_function[0]) );
			//R.push_back( m * (b - x_k) + ( (objective_function[2] - objective_function[1]) * (objective_function[2] - objective_function[1]) ) / (m * (b - x_k)) - 2 * (objective_function[2] + objective_function[1]) );
			//R.push_back(f_R(m, a, x_k, objective_function[0], objective_function[1]));
			//R.push_back(f_R(m, x_k, b, objective_function[1], objective_function[2]));