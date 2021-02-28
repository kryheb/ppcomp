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
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

#include "gnuplot_i.h"

#define RAD(deg) (deg) * M_PI / 180
#define DEG(rad) (rad) * 180 / M_PI

using namespace std::literals;

struct Point {
  double x{};
  double y{};
  double z{};
};

std::ostream& operator<<(std::ostream& os, const Point& p) {
  os << p.x << " " << p.y << " " << p.z << '\n';
  return os;
}

struct PositionProperties {
  double N{};
  double i{};
  double w{};
  double a{};
  double e{};
  double M{};
  double oblec{};
};

struct Position {
  std::string gnuplotFormula(int day) {
    auto filename = "datafile_"s + label() + ".dat";
    {
      std::ofstream f;
      f.open(filename, std::ios::out | std::ios::trunc);
      f << heliocentricPosition(day);
      f.close();
    }
    std::stringstream sstr;
    sstr << "\"" << filename << "\" using 1:2:3:(sprintf(\"" << label() << "\"))"
         << " with labels point lc " << color() << " pointtype 7 pointsize 1 offset char 1,1 notitle";
    return sstr.str();
  }

protected:
  virtual int color() const = 0;
  virtual std::string label() const = 0;
  virtual const PositionProperties properties(int day) const = 0;

private:
  std::string xyz_position(int day) {
    auto pp = properties(day);
    double L{[&] {
      auto result = std::fmod(pp.w + pp.M, 360);
      return result < 0 ? result + 360 : result;
    }()};
    double E{pp.M + DEG(pp.e) * sin(pp.M * M_PI / 180) * (1 + pp.e * cos(pp.M * M_PI / 180))};
    double x{pp.a * (cos(E * M_PI / 180) - pp.e)}, y{pp.a * sin(E * M_PI / 180) * sqrt(1 - pp.e * pp.e)};
    double r{sqrt(x * x + y * y)};
    double v{atan2(y, x) * 180 / M_PI};
    double lon{v + pp.w};
    double _x{r * cos(lon * M_PI / 180)};
    double _y{r * sin(lon * M_PI / 180)};
    double _z{0};
    double xequat{_x};
    double yequat{_y * cos(pp.oblec * M_PI / 180) - 0.0 * sin(pp.oblec * M_PI / 180)};
    double zequat{_y * sin(pp.oblec * M_PI / 180) - 0.0 * cos(pp.oblec * M_PI / 180)};
    double _r{sqrt(xequat * xequat + yequat * yequat + zequat * zequat)};
    double RA{atan2(yequat, xequat) * 180 / M_PI};
    double Decl{atan2(zequat, sqrt(xequat * xequat + yequat * yequat)) * 180 / M_PI};

    std::stringstream sstr;
    sstr << xequat << " " << yequat << " " << zequat;
    return sstr.str();
  }

  double approx_E(const PositionProperties& pp) {
    double E0 = pp.M + pp.e * sin(RAD(pp.M)) * (1 + pp.e * cos(RAD(pp.M)));
    double E1 = E0 - (E0 - pp.e * sin(RAD(E0)) - pp.M) / (1 - pp.e * cos(RAD(E0)));
    while (E0 - E1 >= 0.005) {
      E0 = E1;
      E1 = E0 - (E0 - pp.e * sin(RAD(E0)) - pp.M) / (1 - pp.e * cos(RAD(E0)));
    }
    return E1;
  }

  Point heliocentricPosition(int day) {
    auto pp = properties(day);

    double E{approx_E(pp)};
    double x{pp.a * (cos(RAD(E)) - pp.e)};
    double y{pp.a * sin(RAD(E)) * sqrt(1 - pp.e * pp.e)};
    double r{sqrt(x * x + y * y)};
    double v{DEG(atan2(y, x))};

    double xeclip{r * (cos(RAD(pp.N)) * cos(RAD(v + pp.w)) - sin(RAD(pp.N)) * sin(RAD(v + pp.w)) * cos(RAD(pp.i)))};
    double yeclip{r * (sin(RAD(pp.N)) * cos(RAD(v + pp.w)) + cos(RAD(pp.N)) * sin(RAD(v + pp.w)) * cos(RAD(pp.i)))};
    double zeclip{r * sin(RAD(v + pp.w)) * sin(RAD(pp.i))};

    return Point{xeclip, yeclip, zeclip};
  }
};

struct Sun : Position {
  virtual int color() const { return 1; }
  std::string label() const override { return "Sun"; }
  const PositionProperties properties(int day) const override {
    return {.N{0},
            .i{0},
            .w{282.9404 + 4.70935E-5 * day},
            .a{0},
            .e{0.016709 + 1.151E-9 * day},
            .M{[day] {
              auto result = std::fmod(356.0470 + 0.9856002585 * day, 360);
              return result < 0 ? result + 360 : result;
            }()},
            .oblec{23.4393 - 3.563E-7 * day}};
  }
};

struct Mercury : Position {
  virtual int color() const { return 2; }
  std::string label() const override { return "Mercury"; }
  const PositionProperties properties(int day) const override {
    return {.N{48.3313 + 3.24587E-5 * day},
            .i{7.0047 + 5.00E-8 * day},
            .w{29.1241 + 1.01444E-5 * day},
            .a{0.387098},
            .e{0.205635 + 5.59E-10 * day},
            .M{[day] {
              auto result = std::fmod(168.6562 + 4.0923344368 * day, 360);
              return result < 0 ? result + 360 : result;
            }()},
            .oblec{23.4393 - 3.563E-7 * day}};
  }
};

struct Mars : Position {
  virtual int color() const { return 3; }
  std::string label() const override { return "Mars"; }
  const PositionProperties properties(int day) const override {
    return {.N{49.5574 + 2.46590E-5 * day},
            .i{1.8497 - 1.78E-8 * day},
            .w{286.5016 + 2.92961E-5 * day},
            .a{1.523688},
            .e{0.093405 + 2.526E-9 * day},
            .M{[day] {
              auto result = std::fmod(18.6021 + 0.5240207766 * day, 360);
              return result < 0 ? result + 360 : result;
            }()},
            .oblec{23.4393 - 3.563E-7 * day}};
  }
};

struct Jupiter : Position {
  virtual int color() const { return 4; }
  std::string label() const override { return "Jupiter"; }
  const PositionProperties properties(int day) const override {
    return {.N{100.4542 + 2.76854E-5 * day},
            .i{1.3030 + 1.557E-7 * day},
            .w{273.8777 + 1.64505E-5 * day},
            .a{5.20256},
            .e{0.048498 + 4.469E-9 * day},
            .M{[day] {
              auto result = std::fmod(19.8950 + 0.0830853001 * day, 360);
              return result < 0 ? result + 360 : result;
            }()},
            .oblec{23.4393 - 3.563E-7 * day}};
  }
};

int main() {
  gnuplot_ctrl* h1;

  h1 = gnuplot_init();
  gnuplot_cmd(h1, "set xrange [-10:10]");
  gnuplot_cmd(h1, "set yrange [-10:10]");
  gnuplot_cmd(h1, "set zrange [-5:5]");
  gnuplot_cmd(h1, "set ticslevel 0");

  Sun sun{};
  Mercury mercury{};
  Mars mars{};
  Jupiter jupiter{};

  std::vector<std::unique_ptr<Position>> positions{};
  positions.emplace_back(std::make_unique<Sun>());
  positions.emplace_back(std::make_unique<Mercury>());
  positions.emplace_back(std::make_unique<Mars>());
  positions.emplace_back(std::make_unique<Jupiter>());

  for (int i = 0;; ++i) {
    std::stringstream sstr;
    sstr << "splot ";
    for (const auto& p : positions) {
      sstr << p->gnuplotFormula(i) << ", ";
    }

    gnuplot_cmd(h1, sstr.str().c_str());

    usleep(10000);
  }

  for (;;)
    ;

  gnuplot_close(h1);

  return 0;
}
