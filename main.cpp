/*
 *  Copyright (C) 2021  Krystian Heberlein <krystianheberlein@gmail.com>
 *
 *  This program visualizes the planetary positions in the solar system on the given day.
 *  The plot visualization using Gnuplot http://www.gnuplot.info/
 *  and Gnuplot C interface http://ndevilla.free.fr/gnuplot/.
 *  The computations were made based on the equations from
 *  http://www.stjarnhimlen.se/comp/ppcomp.html
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <unistd.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>
#include <unordered_map>

#include "gnuplot_i.h"

#define RAD(deg) (deg) * M_PI / 180
#define DEG(rad) (rad) * 180 / M_PI

using namespace std::literals;

auto constexpr daysInCentury{36525};
struct OrbitalElement {
  double value{};
  double centennialRate{}; 

  inline double valueOnDay(int day) const { return value + (double)(centennialRate  / daysInCentury) * (double)day; }
};

struct OrbitalElements {
  OrbitalElement a{};
  OrbitalElement e{};
  OrbitalElement i{};
  OrbitalElement Omega{};
  OrbitalElement omega{};
  OrbitalElement L{};
};


static const std::unordered_map<std::string, OrbitalElements> orbitalElements{
/*
 * +--------------+-------------+-------------+------------+-------------+---------+--------------+------------+--------------+------------+---------------+-----------+--------------------+
 * | planet       |     a[AU]   |    a/cy     |     e      |    e/cy     | i[deg]  |     i"/cy    | Omega[deg] |   Omega"/cy  | omega[deg] |    omega"/cy  |  L [deg]  |       L"/cy        | 
 * +--------------+-------------+-------------+------------+-------------+---------+--------------+------------+--------------+------------+---------------+-----------+--------------------+
 */ 
    {"Mercury",   {{ 0.38709893,  0.00000066}, {0.20563069,  0.00002527}, { 7.00487, -23.51/3600}, { 48.33167,   -446.30/3600}, { 77.45645,   573.57/3600}, {252.25084, 538101628.29/3600}}},
    {"Venus",     {{ 0.72333199,  0.00000092}, {0.00677323, -0.00004938}, { 3.39471,  -2.86/3600}, { 76.68069,   -996.89/3600}, {131.53298,  -108.80/3600}, {181.97973, 210664136.06/3600}}},
    {"Earth",     {{ 1.00000011, -0.00000005}, {0.01671022, -0.00003804}, { 0.00005, -46.94/3600}, {-11.26064, -18228.25/3600}, {102.94719,  1198.28/3600}, {100.46435, 129597740.63/3600}}},
    {"Mars",      {{ 1.52366231, -0.00007221}, {0.09341233,  0.00011902}, { 1.85061, -25.47/3600}, { 49.57854,  -1020.19/3600}, {336.04084,  1560.78/3600}, {355.45332,  68905103.78/3600}}},
    {"Jupiter",   {{ 5.20336301,  0.00060737}, {0.04839266, -0.00012880}, { 1.30530,  -4.15/3600}, {100.55615,   1217.17/3600}, { 14.75385,   839.93/3600}, { 34.40438,  10925078.35/3600}}},
    {"Saturn",    {{ 9.53707032, -0.00301530}, {0.05415060, -0.00036762}, { 2.48446,   6.11/3600}, {113.71504,  -1591.05/3600}, { 92.43194, -1948.89/3600}, { 49.94432,   4401052.95/3600}}},
    {"Uranus",    {{19.19126393,  0.00152025}, {0.04716771, -0.00019150}, { 0.76986,  -2.09/3600}, { 74.22988,  -1681.40/3600}, {170.96424,  1312.56/3600}, {313.23218,   1542547.79/3600}}},
    {"Neptune",   {{30.06896348, -0.00125196}, {0.00858587,  0.0000251 }, { 1.76917,  -3.64/3600}, { 31.72169,   -151.25/3600}, { 44.97135,  -844.43/3600}, {304.88003,    786449.21/3600}}}
};

double approxE(double M, double e) {
  double E0 = M + e * sin(RAD(M)) * (1 + e * cos(RAD(M)));
  double E1 = E0 - (E0 - e * sin(RAD(E0)) - M) / (1 - e * cos(RAD(E0)));
  while (E0 - E1 >= 0.005) {
    E0 = E1;
    E1 = E0 - (E0 - e * sin(RAD(E0)) - M) / (1 - e * cos(RAD(E0)));
  }
  return E1;
}

struct Point {
  double x{};
  double y{};
  double z{};
};

Point heliocentricPosition(const OrbitalElements& oe, int day) {
  double L{oe.L.valueOnDay(day)};
  double omega{oe.omega.valueOnDay(day)};
  // double M{L-omega};
  double M{std::fmod(L-omega, 360)};
  M = M < 0 ? M + 360 : M;
  double e{oe.e.valueOnDay(day)};
  double E{approxE(M, e)};
  double a{oe.a.valueOnDay(day)};
  double Omega{oe.Omega.valueOnDay(day)};
  double i{oe.i.valueOnDay(day)};

  double x{a * (cos(RAD(E)) - e)};
  double y{a * sin(RAD(E)) * sqrt(1 - e * e)};
  double r{sqrt(x * x + y * y)};
  double v{DEG(atan2(y, x))};

  double xeclip{r * (cos(RAD(Omega)) * cos(RAD(v + omega)) - sin(RAD(Omega)) * sin(RAD(v + omega)) * cos(RAD(i)))};
  double yeclip{r * (sin(RAD(Omega)) * cos(RAD(v + omega)) + cos(RAD(Omega)) * sin(RAD(v + omega)) * cos(RAD(i)))};
  double zeclip{r * sin(RAD(v + omega)) * sin(RAD(i))};

  return Point{xeclip, yeclip, zeclip};
}


std::ostream& operator<<(std::ostream& os, const Point& p) {
  os << p.x << " " << p.y << " " << p.z << '\n';
  return os;
}

struct GnuplotProperties {
  std::string color{};
  int size{};
};

static const std::unordered_map<std::string, GnuplotProperties> gnuplotProperties{
  {"Mercury",   {"#504e51", 1}},
  {"Venus",     {"#eed053", 2}},
  {"Earth",     {"green",   2}},
  {"Mars",      {"#bc2731", 2}},
  {"Jupiter",   {"#e36e4b", 4}},
  {"Saturn",    {"#ab604a", 4}},
  {"Uranus",    {"#4fD0e7", 4}},
  {"Neptune",   {"#4b70dd", 4}}
};

std::string gnuplotFormula(std::string planet, int day) {
  auto positionFileName = "./ppcomp_data/position_"s + planet + ".dat";
  auto orbitFileName = "./ppcomp_data/orbit_"s + planet + ".dat";
  auto position = heliocentricPosition(orbitalElements.at(planet), day);
  {
    std::ofstream f;
    f.open(positionFileName, std::ios::out | std::ios::trunc);
    f << position;
    f.close();
  }
  {
    std::ofstream f;
    f.open(orbitFileName, std::ios_base::app);
    f << position;
    f.close();
  }
  auto properties = gnuplotProperties.at(planet);
  std::stringstream sstr;
  sstr << "\"" << positionFileName << "\" using 1:2:3:(sprintf(\"" << planet << "\"))"
        << " with labels point lc rgb '" << properties.color << "' pointtype 7 pointsize " << properties.size
        << " offset char 1,1 notitle,"
        << "\"" << orbitFileName << "\" using 1:2:3 "
        << "with lines lc \"grey\"  notitle";
  return sstr.str();
}


int main() {
  gnuplot_ctrl* h1;

  h1 = gnuplot_init();
  gnuplot_cmd(h1, "set xrange [-10:10]");
  gnuplot_cmd(h1, "set yrange [-10:10]");
  gnuplot_cmd(h1, "set zrange [-5:5]");
  gnuplot_cmd(h1, "set ticslevel 0");

  std::filesystem::remove_all("./ppcomp_data");
  std::filesystem::create_directory("./ppcomp_data");

  for (int i = 0; i < 10000; ++i) {

    std::stringstream sstr;
    sstr << "splot ";
    for (const auto& planet : orbitalElements) {
       sstr << gnuplotFormula(planet.first, i) << ", ";
    }

    gnuplot_cmd(h1, sstr.str().c_str());

    usleep(50000);
  }

  gnuplot_close(h1);


  return 0;
}

