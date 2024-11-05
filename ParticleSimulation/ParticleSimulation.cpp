#include <iostream>
#include "RK4Stepper.h"

void harmonicOscillator(std::array<double, 2>& x, std::array<double, 2>& dxdt, double t) {
	dxdt[0] = x[1];
	dxdt[1] = -x[0];
}

int main()
{
	std::array<double, 2> x0 = { 3, 0 };
	RK4Stepper<double, 2> stepper(harmonicOscillator, x0);

	double stopT = 151 * 3.141;
	double dt = 0.0001;
	int stepsCnt = (int)(stopT / dt);
	for (int i = 0; i < stepsCnt; i++) {
		stepper.doStep(dt);
	}
	
	auto& x = stepper.getX();
	std::cout << x[0] << ' ' << x[1] << '\n';
}
