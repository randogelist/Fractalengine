#pragma once
#include "quantum_math.hpp"
#include <cstdint>
#include <limits>
#include <ostream>
#include <string>

// Classical dissipative kicked top on the UNIT sphere. One kick is
// F = D_tau o T_k o R_y(beta), with T_k(q)=R_z(k*q.z)q.
// D_tau is the time-tau superradiant damping flow on S^2, not radial
// contraction into the Bloch ball. Every finite-tau map is invertible.
// References: Carlo et al., arXiv:1809.01716, Eqs. (4)-(6);
// Braun, arXiv:chao-dyn/9910009, Fig. 1, (k,beta,tau)=(4,2,.5).
// The latter is a published numerical attractor example, not a theorem
// proving fractality or a particular Hausdorff dimension.
namespace spin {
using Vec3 = quantum::Vec3;

struct Settings {
    double twist=4, tilt=2, damping=.5;
    int burnPower=10, windowPower=8, gridPower=7;
    double kMin=0, kMax=12, tauMin=0, tauMax=1.5;
};

inline bool valid(const Settings& s) {
    auto range=[](double x,double a,double b) {
        return std::isfinite(x) && x>=a && x<=b;
    };
    return range(s.twist,0,20) && range(s.tilt,0,quantum::pi) &&
           range(s.damping,0,2) && s.burnPower>=6 && s.burnPower<=13 &&
           s.windowPower>=5 && s.windowPower<=11 &&
           s.gridPower>=5 && s.gridPower<=8 &&
           range(s.kMin,0,20) && range(s.kMax,0,20) &&
           range(s.tauMin,0,2) && range(s.tauMax,0,2) &&
           s.kMax-s.kMin>=1e-6 && s.tauMax-s.tauMin>=1e-6;
}

// Deterministic unit tangent, including at the north/south poles.
inline Vec3 tangentAt(Vec3 r) {
    const Vec3 axis=std::abs(r.z)<.8 ? Vec3{0,0,1} : Vec3{0,1,0};
    return quantum::normalized(quantum::cross(axis,r));
}

// Only initial states use pseudorandom numbers. The dynamics itself has
// neither random branch choices nor reseeding. Streams are keyed by chain
// identity, so scheduling and thread count do not select different starts.
inline uint64_t mix64(uint64_t& state) {
    uint64_t z=(state+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}
inline Vec3 initialPoint(uint32_t seed,int chainId) {
    uint64_t state=(uint64_t(seed)<<32)^uint64_t(uint32_t(chainId))^
                   UINT64_C(0x8b8b8b8b52eaf739);
    const double a=double(mix64(state)>>11)*(1.0/9007199254740992.0);
    const double b=double(mix64(state)>>11)*(1.0/9007199254740992.0);
    const double z=2*a-1, phi=2*quantum::pi*b;
    const double radius=std::sqrt(std::max(0.0,1-z*z));
    return quantum::normalized({radius*std::cos(phi),radius*std::sin(phi),z});
}

struct Step {
    Vec3 point, tangent;
    double logStretch=0, logArea=0;
};

// Precomputed constants are shared by all steps of a fixed-parameter orbit.
// Inputs to map/step are unit points; step's tangent is projected and
// normalized before propagation, and its output tangent is unit length.
class Model {
    double k_, c_, s_, t_, sech_, logCosh_;

    Vec3 rotateY(Vec3 r) const {
        return {c_*r.x+s_*r.z,r.y,-s_*r.x+c_*r.z};
    }
    static Vec3 rotateZ(Vec3 r,double c,double s) {
        return {c*r.x-s*r.y,s*r.x+c*r.y,r.z};
    }
    Vec3 damp(Vec3 u,double denominator) const {
        return {sech_*u.x/denominator,sech_*u.y/denominator,
                (u.z-t_)/denominator};
    }

public:
    explicit Model(const Settings& settings) {
        if (!valid(settings)) throw std::invalid_argument("Invalid spin settings");
        k_=settings.twist;
        c_=std::cos(settings.tilt); s_=std::sin(settings.tilt);
        t_=std::tanh(settings.damping);
        const double coshTau=std::cosh(settings.damping);
        sech_=1/coshTau; logCosh_=std::log(coshTau);
    }

    Vec3 map(Vec3 r) const {
        const Vec3 q=rotateY(r);
        const double angle=k_*q.z, c=std::cos(angle), s=std::sin(angle);
        const Vec3 u=rotateZ(q,c,s);
        const double denominator=1-t_*u.z;
        // For unit inputs and 0<=tau<=2 this denominator is >=1-tanh(2)>0.
        // No clipping, bailout, pole special case or artificial reinjection.
        return quantum::normalized(damp(u,denominator));
    }

    Step step(Vec3 r,Vec3 tangent) const {
        tangent=quantum::normalized(tangent-r*quantum::dot(r,tangent));
        const Vec3 q=rotateY(r), v=rotateY(tangent);
        const double angle=k_*q.z, c=std::cos(angle), s=std::sin(angle);
        const Vec3 u=rotateZ(q,c,s);
        // d T_k(q)[v] = R_z(k*q.z) [v + k*v.z*(e_z cross q)].
        const Vec3 w=rotateZ(v+Vec3{-q.y,q.x,0}*(k_*v.z),c,s);
        const double denominator=1-t_*u.z;
        const double factor=sech_/denominator;
        const Vec3 rawPoint=damp(u,denominator);
        const double pointNorm=quantum::length(rawPoint);
        const Vec3 point=rawPoint/pointNorm;
        // Ambient differential of D_tau, restricted to the input tangent:
        // dD_xy = (sech/den)*(w_xy + u_xy*tanh*w_z/den)
        // dD_z  = sech^2*w_z/den^2.
        const double correction=t_*w.z/denominator;
        Vec3 derivative={factor*(w.x+u.x*correction),
                         factor*(w.y+u.y*correction),factor*factor*w.z};
        // Account for the final roundoff normalization of the point. In
        // exact arithmetic pointNorm=1 and derivative is already tangent.
        derivative=(derivative-point*quantum::dot(point,derivative))/pointNorm;
        const double stretch=quantum::length(derivative);
        if (!(stretch>0) || !std::isfinite(stretch) || !quantum::finite(point))
            throw std::runtime_error("Invalid spin tangent propagation");
        // R_y and the latitude twist preserve spherical area. Damping is
        // conformal with scale sech(tau)/(1-u.z*tanh(tau)); its area
        // factor can EXCEED one near the north pole despite its name.
        return {point,derivative/stretch,std::log(stretch),
                -2*logCosh_-2*std::log(denominator)};
    }
};

inline Vec3 map(Vec3 r,const Settings& s) { return Model(s).map(r); }
inline Step step(Vec3 r,Vec3 tangent,const Settings& s) {
    return Model(s).step(r,tangent);
}

// Numerical contracts supplied as SOURCE for the user's native selftest.
// They have not been run in the package-authoring environment. Passing
// contracts checks floating-point implementation, not attractor existence,
// asymptotic Lyapunov signs, topology, convergence or dimension.
inline int selfTest(std::ostream& out) {
    int checks=0;
    auto require=[&](bool condition,const char* name) {
        ++checks;
        if (!condition) {
            out << "SPIN_CHECK_FAIL=" << name << '\n';
            throw std::runtime_error(std::string("Spin self-test: ")+name);
        }
    };
    Settings s;
    require(valid(s),"settings_default");
    for (double k : {0.0,20.0}) for (double beta : {0.0,quantum::pi})
        for (double tau : {0.0,2.0}) {
            Settings edge=s; edge.twist=k; edge.tilt=beta; edge.damping=tau;
            require(valid(edge),"settings_parameter_edges");
        }
    Settings edge=s;
    edge.burnPower=6; edge.windowPower=5; edge.gridPower=5;
    require(valid(edge),"settings_budget_low");
    edge.burnPower=13; edge.windowPower=11; edge.gridPower=8;
    edge.kMax=20; edge.tauMax=2;
    require(valid(edge),"settings_budget_high");
    Settings bad=s; bad.twist=20.01;
    require(!valid(bad),"reject_twist_outside_range");
    bad=s; bad.damping=-.001;
    require(!valid(bad),"reject_negative_damping");
    bad=s; bad.tilt=std::numeric_limits<double>::quiet_NaN();
    require(!valid(bad),"reject_nonfinite_tilt");
    bad=s; bad.kMin=bad.kMax;
    require(!valid(bad),"reject_degenerate_k_interval");
    bad=s; bad.tauMin=bad.tauMax;
    require(!valid(bad),"reject_degenerate_tau_interval");
    bad=s; bad.gridPower=9;
    require(!valid(bad),"reject_grid_budget");
    const std::vector<Vec3> points={{1,0,0},{0,1,0},{0,0,1},{0,0,-1},
        quantum::normalized({.31,-.52,.79}),initialPoint(1,0),initialPoint(9,13)};
    require(quantum::length(initialPoint(1,0)-initialPoint(1,0))==0,
            "seed_repeatability");
    for (double k : {0.0,4.0,20.0}) for (double beta : {0.0,2.0,quantum::pi})
        for (double tau : {0.0,.5,2.0}) {
            s.twist=k; s.tilt=beta; s.damping=tau;
            const Model model(s);
            for (const Vec3 r : points) {
                const Vec3 v=tangentAt(r), v2=quantum::cross(r,v);
                const Step a=model.step(r,v), b=model.step(r,v2);
                require(quantum::finite(a.point)&&quantum::finite(a.tangent)&&
                        std::abs(quantum::length(a.point)-1)<3e-13&&
                        std::abs(quantum::length(a.tangent)-1)<3e-13,
                        "unit_point_and_tangent");
                require(std::abs(quantum::dot(a.point,a.tangent))<3e-13,
                        "tangent_orthogonal");
                require(quantum::length(model.map(r)-a.point)<3e-13,
                        "map_step_consistency");
                const double h=1e-6;
                const Vec3 plus=model.map(r*std::cos(h)+v*std::sin(h));
                const Vec3 minus=model.map(r*std::cos(h)-v*std::sin(h));
                const Vec3 numeric=(plus-minus)/(2*h);
                const Vec3 analytic=a.tangent*std::exp(a.logStretch);
                require(quantum::length(numeric-analytic)<
                        3e-5*(1+quantum::length(analytic)),
                        "analytic_tangent_vs_geodesic_difference");
                const double determinant=quantum::dot(a.point,
                    quantum::cross(a.tangent,b.tangent))*
                    std::exp(a.logStretch+b.logStretch);
                require(determinant>0&&std::abs(std::log(determinant)-a.logArea)<2e-9,
                        "oriented_spherical_area_jacobian");
                if (tau==0) require(std::abs(a.logArea)<3e-13,
                                    "undamped_area_preservation");
                if (k==0) require(std::abs(2*a.logStretch-a.logArea)<2e-11,
                                  "zero_twist_conformal");
                if (k==0&&tau==0) {
                    const Vec3 rotation={std::cos(beta)*r.x+std::sin(beta)*r.z,
                        r.y,-std::sin(beta)*r.x+std::cos(beta)*r.z};
                    require(quantum::length(a.point-rotation)<3e-13&&
                            std::abs(a.logStretch)<3e-13,"rotation_limit");
                }
            }
        }
    s=Settings{}; s.twist=0; s.tilt=0;
    const Model pureDamping(s);
    require(quantum::length(pureDamping.map({0,0,1})-Vec3{0,0,1})<3e-13&&
            quantum::length(pureDamping.map({0,0,-1})-Vec3{0,0,-1})<3e-13,
            "damping_poles_fixed");
    require(pureDamping.step({0,0,1},{1,0,0}).logArea>0&&
            pureDamping.step({0,0,-1},{1,0,0}).logArea<0,
            "north_expands_south_contracts");
    out << "SPIN_NUMERICAL_CONTRACTS=PASS\n";
    return checks;
}
} // namespace spin
