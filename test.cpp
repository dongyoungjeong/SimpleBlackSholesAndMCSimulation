#include <iostream>
using namespace std;

const int NUMBER_OF_COUPONS = 50;
const double COUPON_RATE = .1;

int main() {
	double accrued = 1 / pow(1.03, NUMBER_OF_COUPONS);
	double sum = 0;
	for (int i = 1; i <= NUMBER_OF_COUPONS; i++) sum += (COUPON_RATE / 2) / pow(1.03, i);
	cout << sum + accrued << endl;


	return 0;
}