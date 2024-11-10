#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
using namespace std;
#define _USE_MATH_DEFINES

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
		os << "arg_extremum = " << obj.x << "\nobj_func_extremum = " << obj.y
		   << "\ncount operation to search = " << obj.k << endl << endl;
		return os;
	}
};

class Segment {
public:
	double a;
	double b;

	Segment(): a(), b(){}
	Segment(double a_, double b_) {
		a = a_;
		b = b_;
	}
};

class SegmentForQueue: public Segment{
public:
	double R;
	Segment seg_vec;

	SegmentForQueue(double a_, double b_, double R_, Segment& seg_vec_) {
		seg_vec = seg_vec_;
		a = a_;
		b = b_;
		R = R_;	
	}
	bool operator > (const SegmentForQueue& obj) {
		return R > obj.R;
	}
	bool operator < (const SegmentForQueue& obj) {
		return R < obj.R;
	}
	
};

struct Comp {
	bool operator () (const SegmentForQueue& obj1, const SegmentForQueue& obj2) {
		return obj1.R < obj2.R;
	}
};

Result algorithm_for_searching(double a, double b, double (*f)(double x), double r, int max_count_operation, double eps) {
	int count_operation = 0;
	vector<Segment>  seg_v;
	priority_queue<SegmentForQueue, vector<SegmentForQueue>, Comp> seg_q;
	double res_arg;
	if (f(a) > f(b)) {
		res_arg = b;
	}
	else {
		res_arg = a;
	}

	double max_M = abs((f(b) - f(a)) / (b - a)), m = 0;
	if (max_M == 0)
		m = 1;
	else
		m = r * max_M;
	

	double x_k = 0.5 * (b + a) - (f(b) - f(a)) / (2 * m);
	if (f(res_arg) > f(x_k)) {
		res_arg = x_k;
	}
	seg_v.push_back(Segment(a, x_k));
	seg_v.push_back(Segment(x_k, b));
	double M1 = abs((f(seg_v[0].b) - f(seg_v[0].a)) / (seg_v[0].b - seg_v[0].a));
	double M2 = abs((f(seg_v[1].b) - f(seg_v[1].a)) / (seg_v[1].b - seg_v[1].a));
	if (M1 > m || M2 > m) {
		m = max(M1, M2);
	}
	double R = f_R(m, seg_v[0].a, seg_v[0].b, f(seg_v[0].a), f(seg_v[0].b));
	seg_q.push(SegmentForQueue(seg_v[0].a, seg_v[0].b, R, seg_v[0]));

	R = f_R(m, seg_v[1].a, seg_v[1].b, f(seg_v[1].a), f(seg_v[1].b));
	seg_q.push(SegmentForQueue(seg_v[1].a, seg_v[1].b, R, seg_v[1]));

	while (count_operation < max_count_operation) {
		SegmentForQueue tmp = seg_q.top();
		seg_q.pop();
		double x_k = 0.5 * (tmp.b + tmp.a) - (f(tmp.b) - f(tmp.a)) / (2 * m);
		if (f(res_arg) > f(x_k)) {
			res_arg = x_k;
		}
		M1 = abs((f(x_k) - f(tmp.a)) / (x_k - tmp.a));
		M2 = abs((f(tmp.b) - f(x_k)) / (f(tmp.b) - x_k));
		tmp.seg_vec = Segment(tmp.a, x_k);
		seg_v.push_back(Segment(x_k, tmp.b));
		if (M1 > m || M2 > m) {
			m = max(M1, M2);
			int s = seg_q.size();
			for (int i = 0; i < s; i++) {
				seg_q.pop();
			}

			for (int i = 1; i < seg_v.size(); i++) {
				R = f_R(m, seg_v[i].a, seg_v[i].b, f(seg_v[i].a), f(seg_v[i].b));
				seg_q.push(SegmentForQueue(seg_v[i].a, seg_v[i].b, R, seg_v[i]));
			}
		}
		else {
			R = f_R(m, tmp.a, x_k, f(tmp.a), f(x_k));
			SegmentForQueue new_1 = SegmentForQueue(tmp.a, x_k, R, tmp.seg_vec);
			seg_q.push(new_1);
			R = f_R(m, x_k, tmp.b, f(x_k), f(tmp.b));
			new_1 = SegmentForQueue(x_k, tmp.b, R, seg_v[seg_v.size()-1]);
			seg_q.push(new_1);
		}
		count_operation++;
		if (tmp.b - tmp.a < eps) {
			break;
		}
	}

	Result res = { res_arg, f(res_arg), count_operation};
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
	return -(  (x + sin(x)) * exp(-(x * x))  );
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
	return (x*x - 5*x + 6)/ (x*x +1);
}

double f_8(double x) {
	return (-x + sin(3 * x) - 1);
}

double f_9(double x) {
	return (2 * (x - 3) * (x - 3) + exp(x * x / 2));
}

int main() {
	Result res_f1 = algorithm_for_searching(2.7,  7.5, f_1, 4.29, 100000, 0.00000001);
	Result res_f2 = algorithm_for_searching(0.0, 10.0, f_2, 67,  100000, 0.00000001);
	Result res_f3 = algorithm_for_searching(0.0, 1.2, f_3, 36, 100000, 0.000000001);
	Result res_f4 = algorithm_for_searching(-10.0, 10.0, f_4, 2.5, 100000, 0.000000001);
	Result res_f5 = algorithm_for_searching(2.7, 7.5, f_5, 6, 100000, 0.00000001);
	Result res_f6 = algorithm_for_searching(0.0, 4.0, f_6, 6.5, 100000, 0.000000001); // not true
	Result res_f7 = algorithm_for_searching(-5.0, 5.0, f_7, 6.5, 100000, 0.000000001); // not true
	Result res_f8 = algorithm_for_searching(0.0, 6.5, f_8, 4.0, 100000, 0.000000001);
	Result res_f9 = algorithm_for_searching(-3.0, 3.0, f_9, 85, 100000, 0.000000001);

	cout << res_f1 << endl
		<< res_f2 << endl
		<< res_f3 << endl
		<< res_f4 << endl
		<< res_f5 << endl
		<< res_f6 << endl
		<< res_f7 << endl
		<< res_f8 << endl
		<< res_f9 << endl;
	return 0;
}








