//
// Created by Arian Dovald on 7/9/26.
//

#ifndef PROPAGATE_INTERPOLATE_1D_H
#define PROPAGATE_INTERPOLATE_1D_H

#include <algorithm>
#include <span>
#include <cmath>
#include <vector>

// base class used by all 1D interpolation routines
struct base_interp {
    // integers and pointers
    int n, mm, jsav, cor, dj;
    std::span<const double> xx, yy;
    // constructor
    base_interp(const std::span<const double> x, const std::span<const double> y, const int m)
        : n(x.size()), mm(m), jsav(0), cor(0), xx(x), yy(y) {
        dj = std::min(1, static_cast<int>(std::pow(static_cast<double>(n), 0.25)));
    }
    // given x, return interpolated value from data pointed to by xx and yy
    double interp(const double x) {
        const int jlo = cor ? hunt(x) : locate(x);
        return rawinterp(jlo, x);
    }
    // locate and hunt routines for finding a section of xx data for which x is the midpoint
    int locate(double x);
    int hunt(double x);
    // classes derived from this class must provide the actual interpolation method below
    double virtual rawinterp(int jlo, double x) = 0;
    // destructor
    virtual ~base_interp() = default;
};

// locate routine
inline int base_interp::locate(const double x) {
    if (n < 2 || mm < 2 || mm > n) throw "locate size error during interpolation";
    const bool ascnd = xx[n-1] >= xx[0];
    int jl = 0;
    int ju = n - 1;
    while (ju - jl > 1) {
        const int jm = (ju + jl) >> 1;
        x >= xx[jm] == ascnd ? jl = jm : ju = jm;
    }
    cor = std::abs(jl-jsav) > dj ? 0 : 1;
    jsav = jl;
    return std::max(0, std::min(n-mm, jl-((mm-2)>>1)));
}

// hunt routine
inline int base_interp::hunt(const double x) {
    int jl = jsav, ju;
    if (n < 2 || mm < 2 || mm > n) throw "hunt size error during interpolation";
    const bool ascnd = xx[n-1] >= xx[0];
    if (jl < 0 || jl > n-1) {
        jl = 0;
        ju = n - 1;
    } else {
        int inc = 1;
        if (x >= xx[jl] == ascnd) {
            for (;;) {
                ju = jl + inc;
                if (ju >= n-1) { ju = n-1; break; }
                if ( x < xx[ju] == ascnd) { break; }
                jl = ju; inc += inc;
            }
        } else {
            ju = jl;
            for (;;) {
                ju = jl - inc;
                if (jl <= 0) { jl = 0; break; }
                if (x >= xx[jl] == ascnd) { break; }
                ju = jl; inc += inc;
            }
        }
    }
    while (ju - jl > 1) {
        if (const int jm = (ju + jl) >> 1; x >= xx[jm] == ascnd) {
            jl = jm;
        } else { ju = jm; }
    }
    cor = std::abs(jl-jsav) > dj ? 0 : 1;
    jsav = jl;
    return std::max(0, std::min(n-mm, jl-((mm-2)>>1)));
}

// cubic spline interpolation object
struct spline_interp : base_interp {
    // may hold derivative info
    std::vector<double> y2;
    // constructor
    spline_interp(const std::span<const double> xv, const std::span<const double> yv,
        const double yp1 = 1.e99, const double ypn=1.e99)
        : base_interp(xv, yv,2), y2(xv.size()) { sety2(xv.data(), yv.data(), yp1, ypn); }
    // calculate derivatives routine
    void sety2(const double* xv, const double* yv, double yp1, double ypn);
    // rawinterp routine
    double rawinterp(int jl, double x) override;
};

// second derivatives routine
inline void spline_interp::sety2(const double* xv, const double* yv, const double yp1, const double ypn) {
    double qn, un;
    const int n = y2.size();
    std::vector<double> u(n-1);
    if (yp1 > 0.99e99) { y2[0] = u[0] = 0.0; }
    else { y2[0] = -0.5; u[0] = 3.0/(xv[1]-xv[0])*((yv[1]-yv[0])/(xv[1]-xv[0])-yp1); }
    for (int i = 1; i < n-1; i++) {
        const double sig = (xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]);
        const double p = sig * y2[i-1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1]);
        u[i] = (6.0 * u[i]/(xv[i+1]-xv[i-1]) - sig*u[i-1])/p;
    }
    if (ypn > 0.99e99) { qn = un = 0.0; }
    else { qn = 0.5; un = 3.0/(xv[n-1]-xv[n-2])*(ypn-(yv[n-1]-yv[n-2])/(xv[n-1]-xv[n-2])); }
    y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    for (int k = n-2; k >= 0; k--) { y2[k] = y2[k]*y2[k+1]+u[k]; }
}

// returns cubic spline interpolated value y(x) given pointers to data xx and yy, and stored y2
inline double spline_interp::rawinterp(const int jl, const double x) {
    const int klo = jl, khi = jl + 1;
    const double h = xx[khi] - xx[klo];
    if (h == 0.0) throw "bad input to spline interpolation routine";
    const double a = (xx[khi]-x)/h;
    const double b = (x-xx[klo])/h;
    const double y = a*yy[klo] + b*yy[khi] + ((a*a*a-a)*y2[klo] + (b*b*b-b)*y2[khi])/6.0;
    return y;
}

#endif //PROPAGATE_INTERPOLATE_1D_H
