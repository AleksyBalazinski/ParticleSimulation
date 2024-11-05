#pragma once
#include <functional>
#include <array>

template<typename T, std::size_t N>
class RK4Stepper {
private:
	using SystemFunc = std::function<void(std::array<T, N>&, std::array<T, N>&, double)>;
	SystemFunc system;
	const std::array<T, N>& x0;
	std::array<T, N> x;
	std::array<T, N> dxdt{};
	double currentT;

	std::array<T, N> k1{};
	std::array<T, N> k2{};
	std::array<T, N> k3{};
	std::array<T, N> k4{};

public:
	RK4Stepper(SystemFunc system, std::array<T, N>& x0, double t0 = 0) : system(system), x0(x0) {
		x = x0;
		currentT = t0;
	}

	std::array<T, N>& getX() {
		return x;
	}

	void doStep(double dt) {
		system(x, dxdt, currentT);
		Mul(dxdt, dt, k1);

		system(Add(x, Mul(k1, 0.5, k2), k2), dxdt, currentT + 0.5 * dt);
		Mul(dxdt, dt, k2);

		system(Add(x, Mul(k2, 0.5, k3), k3), dxdt, currentT + 0.5 * dt);
		Mul(dxdt, dt, k3);
	
		system(Add(x, k3, k4), dxdt, currentT + dt);
		Mul(dxdt, dt, k4);

		Mul(k2, 2.0, k2);
		Mul(k3, 2.0, k3);
		Add(Add(Add(k1, k2, k2), k3, k3), k4, k4);
		Mul(k4, 1.0 / 6, k4);
		Add(x, k4, x);
	}

private:
	std::array<T, N>& Mul(std::array<T, N>& a, double s, std::array<T, N>& out) {
		for (int i = 0; i < a.size(); i++) {
			out[i] = s * a[i];
		}

		return out;
	}

	std::array<T, N>& Add(std::array<T, N>& a, std::array<T, N>& b, std::array<T, N>& out) {
		for (int i = 0; i < a.size(); i++) {
			out[i] = a[i] + b[i];
		}

		return out;
	}
};