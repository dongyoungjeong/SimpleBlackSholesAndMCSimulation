#include <iostream>
#include <random>
using namespace std;

static const double kSqrt2 = 1.41421356237309515;
const int POWER = 6;
const int TRIAL = pow(10, POWER);

double normal_cdf(double x, double sigma = 1, double x0 = 0);

int main() {
	double S0 = 42, K = 40, r = .1, sigma = .2, T = .5;
	double d1 = (log(S0 / K) + (r + pow(sigma, 2.) / 2) * T) / sigma / pow(T, .5);
	double d2 = d1 - sigma * pow(T, .5);
	double c = S0 * normal_cdf(d1) - K * exp(-r * T) * normal_cdf(d2);
	double p = K * exp(-r * T) * normal_cdf(-d2) - S0 * normal_cdf(-d1);
//	cout << "c (BSM): " << c << endl;
//	cout << "p (BSM): " << p << endl;

	default_random_engine generator;
	normal_distribution<double> distribution(0., 1.);
	double c_sum = 0, p_sum = 0;
	double c_sum_anti = 0, p_sum_anti = 0;
	double c_mean[POWER] = { 0 }, p_mean[POWER] = { 0 };
	double c_mean_anti[POWER] = { 0 }, p_mean_anti[POWER] = { 0 };
	for (int i = 0; i < TRIAL; i++) {
		double epsilon = distribution(generator);
		double S = S0 * exp((r - pow(sigma, 2) / 2) * T + sigma * epsilon * pow(T, .5));
		double S_anti = S0 * exp((r - pow(sigma, 2) / 2) * T - sigma * epsilon * pow(T, .5));
		for (int j = 1; j < POWER+1; j++) {
			if (i == (int)pow(10, j)) {
				c_mean[j - 1] = c_sum / (double)i;
				p_mean[j - 1] = p_sum / (double)i;
			}
			else if (i == ((int)pow(10, j)/2)) {
				c_mean_anti[j - 1] = c_sum_anti / (double)i / 2.;
				p_mean_anti[j - 1] = p_sum_anti / (double)i / 2.;
			}
		}
		c_sum += max(S - K, 0.) * exp(-r * T);
		p_sum += max(K - S, 0.) * exp(-r * T);
		c_sum_anti += max(S - K, 0.) * exp(-r * T) + max(S_anti - K, 0.) * exp(-r * T);
		p_sum_anti += max(K - S, 0.) * exp(-r * T) + max(K - S_anti, 0.) * exp(-r * T);
	
	}
	c_mean[POWER-1] = c_sum / TRIAL, p_mean[POWER-1] = p_sum / TRIAL;
	for (int i = 0; i < POWER; i++) cout << "c (MC) with " << pow(10, i + 1) << " trials: " << c_mean[i] << "    " << c_mean_anti[i] << endl;
	cout << "c (BSM): " << c << endl;
	cout << endl;
	for (int i = 0; i < POWER; i++) cout << "p (MC) with " << pow(10, i + 1) << " trials: " << p_mean[i] << "    " << p_mean_anti[i] << endl;
	cout << "p (BSM): " << p << endl;

	
	cout << endl << "*** Error Rates ***" << endl;
	for (int i = 0; i < POWER; i++) cout << "c (MC) with " << pow(10, i + 1) << " trials: " << (c-c_mean[i])/c*100. << " %    " << (c-c_mean_anti[i])/c*100. << " %" << endl;
	cout << endl;
	for (int i = 0; i < POWER; i++) cout << "p (MC) with " << pow(10, i + 1) << " trials: " << (p-p_mean[i])/p*100. << " %    " << (p-p_mean_anti[i])/p*100. << " %" <<endl;


	return 0;
}



double normal_cdf(double x, double sigma, double x0) {
	double z = (x - x0) / (sigma * kSqrt2);
	if (z < -1.) return .5 * erfc(-z);
	else return .5 * (1. + erf(z));
}