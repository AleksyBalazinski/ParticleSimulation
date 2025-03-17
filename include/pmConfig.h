#pragma once

enum class InterpolationScheme { NGP, CIC, TSC };
enum class FiniteDiffScheme { TWO_POINT, FOUR_POINT };
enum class GreensFunction { DISCRETE_LAPLACIAN, S1_OPTIMAL, S2_OPTIMAL };