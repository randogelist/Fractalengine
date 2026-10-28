#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

// Conditional weak measurements of a pure qubit, represented on the unit
// Bloch sphere. The finite visit histogram is numerical, while the update and
// channel-mean identities follow directly from the stated Kraus operators.
namespace quantum {
constexpr double pi = 3.1415926535897932384626433832795;
struct Vec3 { double x = 0, y = 0, z = 1; };
inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline Vec3 operator-(Vec3 a) { return {-a.x, -a.y, -a.z}; }
inline Vec3 operator*(Vec3 a, double b) { return {a.x*b, a.y*b, a.z*b}; }
inline Vec3 operator*(double a, Vec3 b) { return b*a; }
inline Vec3 operator/(Vec3 a, double b) { return {a.x/b, a.y/b, a.z/b}; }
inline double dot(Vec3 a, Vec3 b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
inline Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}
inline bool finite(Vec3 v) { return std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z); }
inline double length(Vec3 v) { return std::hypot(std::hypot(v.x,v.y),v.z); }
inline Vec3 normalized(Vec3 v) {
    const double n = length(v);
    if (!(n > 0) || !std::isfinite(n))
        throw std::runtime_error("Invalid quantum vector norm");
    return v/n;
}
struct Settings {
    double strength = .5, kick = 0, axisAzimuth = 0, axisElevation = .6154797086703874;
    double yaw = .4, pitch = .2, zoom = 1, exposure = 1.4;
    int geometry = 0, samplePower = 20;
    uint32_t randomSeed = 1;
};
inline const char* geometryName(int geometry) {
    static const char* names[] = {"Tetrahedron", "Octahedron", "Icosahedron", "Dodecahedron"};
    return names[std::clamp(geometry,0,3)];
}
inline const std::vector<Vec3>& directions(int geometry) {
    static const std::vector<std::vector<Vec3>> all = [] {
        std::vector<std::vector<Vec3>> a(4);
        a[0] = {{1,1,1},{1,-1,-1},{-1,1,-1},{-1,-1,1}};
        a[1] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        const double phi = (1+std::sqrt(5.0))*.5;
        for (int s : {-1,1}) for (int t : {-1,1}) {
            a[2].push_back({0,double(s),t*phi});
            a[2].push_back({double(s),t*phi,0});
            a[2].push_back({s*phi,0,double(t)});
            a[3].push_back({0,s/phi,t*phi});
            a[3].push_back({s/phi,t*phi,0});
            a[3].push_back({s*phi,0,t/phi});
            for (int u : {-1,1}) a[3].push_back({double(s),double(t),double(u)});
        }
        for (auto& geometry : a) for (auto& n : geometry) n = normalized(n);
        return a;
    }();
    return all[std::clamp(geometry,0,3)];
}
inline Vec3 kickAxis(const Settings& s) {
    const double c = std::cos(s.axisElevation);
    return {c*std::cos(s.axisAzimuth), c*std::sin(s.axisAzimuth), std::sin(s.axisElevation)};
}
// Rodrigues rotation. Axis is unit length; a kick is applied BEFORE outcome
// probabilities and measurement. Thus the Kraus family is K_i U, complete.
inline Vec3 rotate(Vec3 r, Vec3 axis, double radians) {
    if (radians == 0) return r;
    const double c=std::cos(radians), s=std::sin(radians);
    return r*c + cross(axis,r)*s + axis*(dot(axis,r)*(1-c));
}
inline double probability(Vec3 r, Vec3 n, double e, int count) {
    const double t=std::clamp(dot(n,r),-1.0,1.0), a=1-e;
    // Equivalent to 1+e^2+2e*t, avoiding cancellation near e=1,t=-1.
    return (a*a + 2*e*(1+t))/(double(count)*(1+e*e));
}
inline Vec3 measure(Vec3 r, Vec3 n, double e) {
    if (e == 0) return r;
    if (e == 1) return n; // The antipodal zero-probability outcome is never sampled.
    const double t=std::clamp(dot(n,r),-1.0,1.0), a=1-e;
    const double den=a*a + 2*e*(1+t);
    // Parallel/perpendicular form of
    // [(1-e^2)r + 2e(1+e*n.r)n]/[1+e^2+2e*n.r].
    const Vec3 q=(r-n*t)*((1-e*e)/den) + n*(1-a*a*(1-t)/den);
    return normalized(q); // Remove accumulated floating-point norm drift.
}
// For every geometry above: sum n_i=0, sum n_i n_i^T=(N/3)I.
// Therefore sum_i K_i^dagger K_i=I and E[r_next|r]=meanFactor(e)*R_kick*r.
inline double meanFactor(double e) { return (1-e*e/3)/(1+e*e); }
inline bool valid(const Settings& s) {
    auto range=[](double v,double lo,double hi) { return std::isfinite(v)&&v>=lo&&v<=hi; };
    return range(s.strength,0,1)&&range(s.kick,-pi,pi)&&range(s.axisAzimuth,-pi,pi)&&
           range(s.axisElevation,-pi*.5,pi*.5)&&range(s.yaw,-pi,pi)&&
           range(s.pitch,-pi*.5,pi*.5)&&range(s.zoom,.4,8)&&range(s.exposure,.25,4)&&
           s.geometry>=0&&s.geometry<=3&&s.samplePower>=14&&s.samplePower<=23;
}
} // namespace quantum
