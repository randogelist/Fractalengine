#pragma once
#include "riemann_witness_math.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// Three explicitly chosen dynamics from the same fixed standardized moments:
// A=1+h+b*b/3, a=1+b*b/4, eta=sqrt(h).
// 1. Anti-polynomial: conj(z)^3-A*conj(z)+c.
// 2. Alternating return: Q_-(Q_+(z)), Q_+/-(z)=z*z-a+/-c.
// 3. Rational pole feedback: z*z+c+eta/z.
// These are new visualization constructions, not deductions of RH.
// All finite critical points are tracked; infinity is always attracting.
// Homogeneous coordinates preserve the genuine pole in generation 3.
// Finite-budget nonescape is unresolved. No connectedness label is inferred.
namespace momenttriple {
using C=std::complex<double>;
using Color=std::array<double,3>;
inline const char *title(int g) {
    return g==1 ? "RH / conjugate Jacobi cubic" : g==2 ?
        "RH / alternating Jacobi return" : "RH / moment pole feedback";
}
inline const char *formula(int g) {
    return g==1 ? "z_next=conj(z)^3-A*conj(z)+c; A=1+h+b^2/3" : g==2 ?
        "z_next=(z^2-a+c)^2-a-c; a=1+b^2/4; one step = two quadratic updates" :
        "z_next=z^2+c+eta/z; eta=sqrt(h)>0; 0 -> infinity -> infinity";
}
enum class Status { Escaped, Unresolved, NumericalLimit, InvalidModel };
inline const char *statusName(Status status) {
    switch(status) {
    case Status::Escaped:return "numerical infinity-basin criterion reached";
    case Status::Unresolved:return "unresolved within iteration budget";
    case Status::NumericalLimit:return "floating-point evaluation limit (not escape)";
    default:return "invalid model or parameter (finite |c| <= 1e6 required)";
    }
}
struct Result {
    Status status=Status::Unresolved;
    C endpoint{},endpointZ{},endpointW{1,0};
    bool atInfinity=false,endpointRepresentable=true;
    int iterations=0;
    double smooth=0,structure=0,phase=0;
};
class Model {
    static bool finite(C z) {return std::isfinite(z.real())&&std::isfinite(z.imag());}
    static bool validParameter(C c) {return finite(c)&&std::abs(c)<=1e6;}
    static void storeEndpoint(Result &r,C Z,C W) {
        r.endpointZ=Z;r.endpointW=W;r.atInfinity=W==C{};
        r.endpointRepresentable=false;
        const double nan=std::numeric_limits<double>::quiet_NaN();
        r.endpoint={nan,nan};
        if(r.atInfinity||!finite(Z)||!finite(W))return;
        using L=std::complex<long double>;
        const L p=L(Z.real(),Z.imag())/L(W.real(),W.imag());
        const long double limit=std::numeric_limits<double>::max();
        if(!std::isfinite(p.real())||!std::isfinite(p.imag())||
           std::abs(p.real())>limit||std::abs(p.imag())>limit)return;
        r.endpoint={static_cast<double>(p.real()),static_cast<double>(p.imag())};
        r.endpointRepresentable=finite(r.endpoint);
    }
    static Color errorColour(Status s) {
        return s==Status::InvalidModel ? Color{{.42,.035,.055}} : Color{{.83,.12,.72}};
    }
    static double cleanHue(double h) {return std::isfinite(h)?h:0;}
    static double cleanBands(double b) {return std::isfinite(b)?std::clamp(b,0.,100.):1;}
public:
    std::shared_ptr<const riemann::Model> source;
    int generation=1,degree=3,infinityPower=3;
    double b=0,h=0,A=0,a=0,eta=0;
    bool validModel=false;
    std::string description,statusText;
    explicit Model(int g):source(riemann::model(6,riemann::Settings{})),generation(g) {
        description=title(g);
        degree=g==2?4:3;infinityPower=g==1?3:g==2?4:2;
        if(g<1||g>3||!source||!source->validModel) {
            statusText="The fixed H2 point-moment source is unavailable";return;
        }
        b=source->standardizedSkew;h=source->standardizedH;
        A=1+h+b*b/3;a=1+b*b/4;eta=std::sqrt(h);
        if(!std::isfinite(b)||!std::isfinite(h)||!(h>0)||
           !std::isfinite(A)||!std::isfinite(a)||!std::isfinite(eta)||!(eta>0)) {
            statusText="Finite standardized b and strictly positive h are required";return;
        }
        validModel=true;
        statusText="Fixed rounded RH moments; finite critical-orbit diagnostics. Nonescape is unresolved.";
    }
    std::vector<C> criticalPoints(C c) const {
        if(generation==1) {const double r=std::sqrt(A/3);return {C(r,0),C(-r,0)};}
        if(generation==2) {const C r=std::sqrt(C(a,0)-c);return {C{},r,-r};}
        const double r=std::cbrt(eta/2);
        return {C(r,0),std::polar(r,2*projects::pi/3),std::polar(r,4*projects::pi/3)};
    }
    C criticalFixedParameter() const {
        if(generation==1) {const double r=std::sqrt(A/3);return C(r+2*A*r/3,0);}
        // Stable form of a+(1-sqrt(1+8*a))/2: makes F_c(0)=0.
        if(generation==2)return C(a-4*a/(1+std::sqrt(1+8*a)),0);
        const double r=std::cbrt(eta/2);return C(r-3*r*r,0);
    }
    double escapeRadius(C c) const {
        if(!validModel||!validParameter(c))return std::numeric_limits<double>::quiet_NaN();
        const double m=std::abs(c);
        const double R=generation==1 ? std::max({2.,std::sqrt(2*(A+1)),std::cbrt(2*m)}) :
            generation==2 ? std::max(2.,std::sqrt(2*(a+m+1))) :
            std::max({2.,2*std::sqrt(m),std::cbrt(4*eta)});
        return (1+1e-12)*R;
    }
    Result orbit(C initial,C c,int budget,std::vector<C>*trail=nullptr) const {
        Result result;
        if(trail)trail->clear();
        result.endpoint=initial;result.endpointRepresentable=finite(initial);
        if(!validModel||!validParameter(c)) {result.status=Status::InvalidModel;return result;}
        if(!finite(initial)) {result.status=Status::NumericalLimit;return result;}
        const double R=escapeRadius(c);
        // Return polynomial: z^4+k2*z^2+k0, with the PLUS phase first.
        const C u=c-a,k2=2.0*u,k0=u*u-a-c;
        const double scale=generation==1 ? std::max({1.,A,std::abs(c)}) :
            generation==2 ? std::max({1.,std::abs(k2),std::abs(k0)}) :
            std::max({1.,std::abs(c),eta});
        if(!std::isfinite(R)||!(R>1)||!std::isfinite(scale)||!finite(k0)) {
            result.status=Status::NumericalLimit;return result;
        }
        const double initialScale=std::max({1.,std::abs(initial.real()),std::abs(initial.imag())});
        C Z=initial/initialScale,W(1/initialScale,0);
        budget=std::clamp(budget,32,1024);
        bool trailOpen=true;
        double sum=0,weights=0,weight=1;
        for(int n=0;n<=budget;++n) {
            result.endpointZ=Z;result.endpointW=W;result.atInfinity=W==C{};
            result.iterations=n;result.smooth=n;
            if(!finite(Z)||!finite(W)||(Z==C{}&&W==C{})) {result.status=Status::NumericalLimit;break;}
            if(trail&&trailOpen) {
                storeEndpoint(result,Z,W);
                if(result.atInfinity||!result.endpointRepresentable||std::abs(result.endpoint)>=1e12)trailOpen=false;
                if(trailOpen)trail->push_back(result.endpoint);
            }
            const double zm=std::abs(Z),wm=std::abs(W),mass=zm*zm+wm*wm;
            const double shape=.55*zm/(zm+wm)+.45*2*std::abs((Z*std::conj(W)).imag())/mass;
            sum+=weight*std::clamp(shape,0.,1.);weights+=weight;weight*=.97;
            if(result.atInfinity||zm/R>wm) {
                result.status=Status::Escaped;
                if(!result.atInfinity) {
                    const double logRadius=std::log(zm)-std::log(wm);
                    result.smooth=std::max(0.,n+1.-std::log(logRadius/std::log(R))/std::log(double(infinityPower)));
                }
                break;
            }
            if(n==budget)break;
            C nextZ,nextW;
            if(generation==1) {
                const C U=std::conj(Z),V=std::conj(W),V2=V*V,V3=V2*V;
                nextZ=(1/scale)*U*U*U-(A/scale)*U*V2+(c/scale)*V3;
                nextW=V3/scale;
            } else if(generation==2) {
                const C Z2=Z*Z,W2=W*W,W4=W2*W2;
                nextZ=(1/scale)*Z2*Z2+(k2/scale)*Z2*W2+(k0/scale)*W4;
                nextW=W4/scale;
            } else {
                const C W2=W*W,product=Z*W2;
                nextZ=(1/scale)*Z*Z*Z+(c/scale)*product+(eta/scale)*W2*W;
                nextW=product/scale;
            }
            const double norm=std::max(std::abs(nextZ),std::abs(nextW));
            if(!finite(nextZ)||!finite(nextW)||!std::isfinite(norm)||!(norm>0)) {
                result.status=Status::NumericalLimit;break;
            }
            const C normalizedZ=nextZ/norm,normalizedW=nextW/norm;
            const bool denominatorMustSurvive=W!=C{}&&(generation!=3||Z!=C{});
            if((denominatorMustSurvive&&normalizedW==C{})||
               (nextZ!=C{}&&normalizedZ==C{})||!finite(normalizedZ)||!finite(normalizedW)) {
                result.status=Status::NumericalLimit;break;
            }
            Z=normalizedZ;W=normalizedW;
        }
        storeEndpoint(result,result.endpointZ,result.endpointW);
        result.structure=weights>0?sum/weights:0;
        if(!result.atInfinity&&result.endpointZ!=C{}&&result.endpointW!=C{}) {
            result.phase=.5+(std::arg(result.endpointZ)-std::arg(result.endpointW))/(2*projects::pi);
            result.phase-=std::floor(result.phase);
        }
        return result;
    }
    Color juliaColour(C initial,C c,int budget,int palette,double hue,double bands) const {
        const auto r=orbit(initial,c,budget);hue=cleanHue(hue);bands=cleanBands(bands);
        if(r.status==Status::InvalidModel||r.status==Status::NumericalLimit)return errorColour(r.status);
        if(r.status==Status::Unresolved)return projects::hsv(hue+.61,.48,.019+.017*r.structure);
        const double stripes=.5+.5*std::cos(bands*(.57*r.smooth+2.4*r.structure));
        const double colour=palette==1 ? hue+.66-.023*r.smooth : palette==2 ?
            hue+.65*r.structure+.18*r.phase : hue+.57+.30*r.phase+.012*r.smooth;
        return projects::hsv(colour,.74,.26+.72*stripes);
    }
    Color atlasColour(C c,int budget,int palette,double hue,double bands) const {
        if(!validModel||!validParameter(c))return errorColour(Status::InvalidModel);
        const auto starts=criticalPoints(c);
        unsigned mask=0;int escaped=0;double smooth=0,shape=0;
        for(size_t k=0;k<starts.size();++k) {
            const auto r=orbit(starts[k],c,budget);
            if(r.status==Status::InvalidModel||r.status==Status::NumericalLimit)return errorColour(r.status);
            shape+=r.structure;
            if(r.status==Status::Escaped) {mask|=1u<<k;++escaped;smooth+=r.smooth;}
        }
        hue=cleanHue(hue);bands=cleanBands(bands);shape/=starts.size();
        if(!escaped)return projects::hsv(hue+.61,.48,.019+.017*shape);
        smooth/=escaped;
        const double stripes=.5+.5*std::cos(bands*(.57*smooth+2.4*shape));
        const double fateHue=.06+.81*mask/double((1u<<starts.size())-1);
        const double colour=palette==1?hue+.66-.023*smooth:palette==2?hue+fateHue+.15*(shape-.5):hue+fateHue;
        return projects::hsv(colour,.75,.25+.73*stripes);
    }
};
inline std::shared_ptr<const Model> model(int generation) {
    static const std::array<std::shared_ptr<const Model>,3> models={{
        std::make_shared<const Model>(1),std::make_shared<const Model>(2),std::make_shared<const Model>(3)}};
    return models[std::clamp(generation,1,3)-1];
}
} // namespace momenttriple
