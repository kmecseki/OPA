#pragma once
#ifndef PHYS_CONSTANTS_H
#define PHYS_CONSTANTS_H

namespace PhysicalConstants {
	constexpr double eps0 = 8.854e-12; // 8.85418782e-12;
	constexpr double hpl = 6.6256e-34; // 6.626068e-34;
	constexpr double c0 = 2.997925000000e8;// 299792458;
	constexpr double hce0 = 0.5 * c0 * eps0;
}

#endif // PHYS_CONSTANTS_H