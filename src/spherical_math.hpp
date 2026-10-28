#pragma once
#include "quantum_math.hpp"
#include <cstddef>
#include <complex>
#include <limits>
#include <ostream>
#include <string>

// Three constructions on S^2, sharing geometry rather than a physical model:
//   1. Four tangent-circle reflections: Apollonian group limit set.
//   2. Four disjoint-circle reflections: Schottky group limit set.
//   3. Uniform backward branches of z^3 + lambda/z^3: Julia measure.
// Finite sampled densities are neither dimension certificates nor, for the
// reflection walks, an assertion of sampling Hausdorff measure.
namespace spherical {
using Vec3 = quantum::Vec3;
using Complex = std::complex<double>;

struct Settings {
    double gap = .025;
    double lambdaLog = -2;
    double lambdaPhase = quantum::pi;
    double boost = 0;
};

inline bool valid(const Settings& s) {
    auto range=[](double x,double a,double b) {
        return std::isfinite(x) && x>=a && x<=b;
    };
    return range(s.gap,0,.38) && range(s.lambdaLog,-8,0) &&
           range(s.lambdaPhase,-quantum::pi,quantum::pi) &&
           range(s.boost,-.85,.85);
}

inline const std::vector<Vec3>& tetrahedron() { return quantum::directions(0); }
inline double tangentOffset() { return 1/std::sqrt(3.0); }
inline double circleOffset(const Settings& s,int generation) {
    return tangentOffset() + (generation==2 ? s.gap : 0.0);
}
inline double reflectionOffset(int generation,const Settings& s) {
    return circleOffset(s,generation);
}

// Reflection in the spherical circle n.r=a, for unit n and |r|=1:
// [(1-a^2)r + 2(a-n.r)n] / [1+a^2-2a*n.r].
// This is an orientation-reversing conformal involution, not a rotation.
// The perpendicular/parallel expression reduces cancellation at circle poles.
inline Vec3 reflect(Vec3 r,Vec3 n,double a) {
    if (!std::isfinite(a) || !(a>=0 && a<1))
        throw std::invalid_argument("Invalid spherical reflection offset");
    const double t=std::clamp(quantum::dot(n,r),-1.0,1.0);
    const double lo=1-a, hi=1+a;
    const double den=lo*lo+2*a*(1-t);
    const double parallel=t>a ? -1+hi*hi*(1-t)/den
                              : 1-lo*lo*(1+t)/den;
    return quantum::normalized((r-n*t)*((1-a*a)/den)+n*parallel);
}
inline Vec3 reflection(Vec3 r,int index,double a) {
    return reflect(r,tetrahedron().at(static_cast<std::size_t>(index)),a);
}
inline Vec3 reflect(Vec3 r,int index,const Settings& s,int generation) {
    return reflection(r,index,circleOffset(s,generation));
}

// A north-axis Mobius coordinate change. It is applied AFTER sampling and
// is not a measurement, kick, or change in branch probabilities. B_-b=B_b^-1.
// Positive strength with a signed axis retains quantum::measure's stable
// positive-strength denominator, including for negative b.
inline Vec3 boost(Vec3 r,double b) {
    if (!std::isfinite(b) || !(b>=-.85 && b<=.85))
        throw std::invalid_argument("Invalid spherical display boost");
    return quantum::measure(r,{0,0,b<0 ? -1.0 : 1.0},std::abs(b));
}
inline Vec3 display(Vec3 r,const Settings& s) { return boost(r,s.boost); }

// Homogeneous coordinates: z=u/v, with [1:0] representing infinity.
// Every routine normalizes its complex pair; there is no pole clipping,
// escape bailout, or stochastic restart in the inverse dynamics.
struct CP1 {
    Complex u{1,0}, v{1,0};
};
inline bool finite(Complex z) {
    return std::isfinite(z.real()) && std::isfinite(z.imag());
}
inline CP1 normalized(CP1 p) {
    const double scale=std::max(std::abs(p.u),std::abs(p.v));
    if (!finite(p.u) || !finite(p.v) || !(scale>0) || !std::isfinite(scale))
        throw std::runtime_error("Invalid spherical homogeneous point");
    p.u/=scale;
    p.v/=scale;
    const double norm=std::hypot(std::abs(p.u),std::abs(p.v));
    p.u/=norm;
    p.v/=norm;
    return p;
}
inline CP1 initialPoint() { return normalized({Complex(.7,.2),Complex(1,0)}); }
inline Complex lambda(const Settings& s) {
    return std::polar(std::pow(10.0,s.lambdaLog),s.lambdaPhase);
}
inline Vec3 sphere(CP1 p) {
    p=normalized(p);
    const Complex c=p.u*std::conj(p.v);
    return quantum::normalized({2*c.real(),2*c.imag(),
                                 std::norm(p.u)-std::norm(p.v)});
}

// The degree-six rational map in homogeneous coordinates:
// [u:v] -> [u^6 + lambda*v^6 : u^3*v^3].
// lambda != 0 throughout the UI range, so these coordinates never vanish
// simultaneously, even at zero or infinity.
inline CP1 forward(CP1 p,const Settings& s) {
    p=normalized(p);
    const Complex u3=p.u*p.u*p.u, v3=p.v*p.v*p.v;
    return normalized({u3*u3+lambda(s)*v3*v3,u3*v3});
}

// For a target [u:v], t=z^3 satisfies v*t^2-u*t+lambda*v=0.
// With q=u +/- sqrt(u^2-4*lambda*v^2), select the sign maximizing |q|.
// The two stable homogeneous roots are [q:2v] and [2lambda*v:q].
// Taking three cube roots of each gives all six inverse branches, counted
// with multiplicity. At a critical value or infinity repeated branches
// must remain repeated: the sampler chooses branch 0..5 uniformly.
inline CP1 inverseMcMullen(CP1 target,int branch,const Settings& s) {
    if (branch<0 || branch>=6)
        throw std::invalid_argument("Invalid McMullen inverse branch");
    target=normalized(target);
    const Complex l=lambda(s);
    const Complex disc=std::sqrt(target.u*target.u-4.0*l*target.v*target.v);
    const Complex plus=target.u+disc, minus=target.u-disc;
    const Complex q=std::abs(plus)>=std::abs(minus) ? plus : minus;
    CP1 t=branch<3 ? CP1{q,2.0*target.v} : CP1{2.0*l*target.v,q};
    t=normalized(t);
    const int k=branch%3;
    const double theta=(std::arg(t.u)+2*quantum::pi*k)/3;
    const double phi=std::arg(t.v)/3;
    return normalized({std::polar(std::cbrt(std::abs(t.u)),theta),
                       std::polar(std::cbrt(std::abs(t.v)),phi)});
}
inline CP1 backward(CP1 target,const Settings& s,int branch) {
    return inverseMcMullen(target,branch,s);
}

// Analytic bound for the actual starting state and allowed lambda range:
// |w_0|<2 and |lambda|<=1. If |z|>=2 then
// |z^3+lambda/z^3| >= 8-1/8 > 2, so all backward iterates stay |z|<2.
// More sharply, each cubic root t obeys |t|<=1+sqrt(2), and the two
// roots have product lambda. Therefore every inverse state has
// |z| >= (10^-8/(1+sqrt(2)))^(1/3) > 0. The homogeneous algorithm additionally
// handles exact zero/infinity inputs, and never needs that bound as a cutoff.

inline double chordalDistance(CP1 a,CP1 b) {
    a=normalized(a);
    b=normalized(b);
    return std::abs(a.u*b.v-a.v*b.u);
}

// Numerical contracts supplied as source for the user's native selftest.
// These checks have not been executed by the package authoring process;
// successful tests establish floating-point contracts, not fractal dimension.
inline int selfTest(std::ostream& out) {
    int checks=0;
    auto require=[&](bool condition,const char* name) {
        ++checks;
        if (!condition) {
            out << "SPHERICAL_CHECK_FAIL=" << name << '\n';
            throw std::runtime_error(std::string("Spherical self-test: ")+name);
        }
    };
    Settings s;
    require(valid(s),"settings_default");
    for (double g : {0.0,.38}) for (double b : {-.85,.85}) {
        Settings edge=s; edge.gap=g; edge.boost=b;
        require(valid(edge),"settings_edges");
    }
    Settings bad=s; bad.gap=-.001;
    require(!valid(bad),"reject_negative_gap");
    bad=s; bad.lambdaLog=-8.001;
    require(!valid(bad),"reject_small_lambda");
    bad=s; bad.lambdaPhase=std::numeric_limits<double>::infinity();
    require(!valid(bad),"reject_nonfinite_phase");
    bad=s; bad.boost=.851;
    require(!valid(bad),"reject_large_boost");
    const auto& ns=tetrahedron();
    for (std::size_t i=0;i<ns.size();++i) for (std::size_t j=i+1;j<ns.size();++j) {
        require(std::abs(quantum::dot(ns[i],ns[j])+1.0/3)<2e-14,"tetrahedral_dot");
        const Vec3 touching=quantum::normalized(ns[i]+ns[j]);
        require(std::abs(quantum::dot(ns[i],touching)-tangentOffset())<2e-14 &&
                std::abs(quantum::dot(ns[j],touching)-tangentOffset())<2e-14,
                "circle_tangency");
        require(quantum::length(reflect(touching,ns[i],tangentOffset())-touching)<2e-13,
                "tangent_point_fixed");
    }
    const std::vector<Vec3> points={{1,0,0},{0,1,0},{0,0,1},{0,0,-1},
        quantum::normalized({.31,-.52,.79}),ns[0],-ns[0]};
    for (double a : {tangentOffset(),tangentOffset()+.025,tangentOffset()+.38}) {
        for (const auto& n : ns) for (const auto& r : points) {
            const Vec3 q=reflect(r,n,a), recovered=reflect(q,n,a);
            require(quantum::finite(q)&&std::abs(quantum::length(q)-1)<2e-13,
                    "reflection_unit_sphere");
            require(quantum::length(recovered-r)<2e-10,"reflection_involution");
        }
    }
    for (double b : {-.85,-.2,0.0,.2,.85}) for (const auto& r : points) {
        const Vec3 q=boost(r,b), recovered=boost(q,-b);
        require(quantum::finite(q)&&std::abs(quantum::length(q)-1)<2e-13,
                "boost_unit_sphere");
        require(quantum::length(recovered-r)<2e-11,"boost_inverse");
    }
    require(quantum::length(sphere({Complex(0,0),Complex(1,0)})-Vec3{0,0,-1})<2e-14,
            "zero_south_pole");
    require(quantum::length(sphere({Complex(1,0),Complex(0,0)})-Vec3{0,0,1})<2e-14,
            "infinity_north_pole");
    for (double log : {-8.0,-2.0,0.0}) for (double phase : {-quantum::pi,-.731,0.0,quantum::pi}) {
        s.lambdaLog=log; s.lambdaPhase=phase;
        const Complex cv=2.0*std::sqrt(lambda(s));
        const std::vector<CP1> targets={initialPoint(),{Complex(0,0),Complex(1,0)},
            {Complex(1,0),Complex(0,0)},{Complex(-.23,.77),Complex(1,0)},
            {cv,Complex(1,0)},{-cv,Complex(1,0)}};
        for (const auto& target : targets) for (int branch=0;branch<6;++branch) {
            const CP1 z=inverseMcMullen(target,branch,s);
            const Vec3 r=sphere(z);
            require(quantum::finite(r)&&std::abs(quantum::length(r)-1)<2e-13,
                    "mcmullen_unit_sphere");
            require(chordalDistance(forward(z,s),target)<2e-10,"mcmullen_inverse_residual");
        }
    }
    out << "SPHERICAL_NUMERICAL_CONTRACTS=PASS\n";
    return checks;
}
} // namespace spherical
