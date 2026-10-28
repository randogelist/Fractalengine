#pragma once
#include "project_math.hpp"
#include "rh_master.hpp"
#include "rh_candidates.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// Native Catalog adapter for finite-Taylor Master / Sub_1 / Sub_2.
// The parameter atlas follows the chosen seed z=0, not all critical points.
// The radius-16 stopping rule is only a display cutoff.  It is NOT a proved
// escape radius for these rational maps.  Numerical failures stay separate.
namespace rhmaster {
using C = std::complex<double>;
using Color = std::array<double, 3>;

constexpr int MapCount = 12;
constexpr int ChartCount = 4;
constexpr int ColouringCount = 8;
struct Experiment {
    int map=0,chart=0,coloring=0,power=2;
    double gain=1.,relaxation=1.,memory=.15,trapRadius=.5,trapAngle=0.,exposure=1.,frequency=1.;
    double twist=0.,blend=.5;
};
inline bool validExperiment(const Experiment& e) {
    auto range=[](double v,double lo,double hi){return std::isfinite(v)&&v>=lo&&v<=hi;};
    return e.map>=0&&e.map<MapCount&&e.chart>=0&&e.chart<ChartCount
        &&e.coloring>=0&&e.coloring<ColouringCount&&e.power>=2&&e.power<=6
        &&range(e.gain,0.,2.)&&range(e.relaxation,.05,1.95)&&range(e.memory,0.,.95)
        &&range(e.trapRadius,.01,3.)&&range(e.trapAngle,-projects::pi,projects::pi)
        &&range(e.exposure,.1,3.)&&range(e.frequency,.2,8.)
        &&range(e.twist,-projects::pi,projects::pi)&&range(e.blend,0.,1.);
}
inline bool sameExperiment(const Experiment& a,const Experiment& b) {
    return a.map==b.map&&a.chart==b.chart&&a.coloring==b.coloring&&a.power==b.power
        &&a.gain==b.gain&&a.relaxation==b.relaxation&&a.memory==b.memory
        &&a.trapRadius==b.trapRadius&&a.trapAngle==b.trapAngle
        &&a.exposure==b.exposure&&a.frequency==b.frequency
        &&a.twist==b.twist&&a.blend==b.blend;
}
inline const char* mapName(int i) {
    static const char* names[MapCount]={"Classic polynomial","Newton + parameter","Conjugate host",
        "Alternating feedback","Inertial feedback","Cayley feedback",
        "Conjugate Cayley","Alternating Cayley","Inertial Cayley","Twisted Cayley",
        "Mixed Tricorn","Mixed Cayley"};
    return i>=0&&i<MapCount?names[i]:"Invalid map";
}
inline const char* chartName(int i) {
    static const char* names[ChartCount]={"Cartesian","Reciprocal","Log-polar","Joukowski"};
    return i>=0&&i<ChartCount?names[i]:"Invalid chart";
}
inline const char* colouringName(int i) {
    static const char* names[ColouringCount]={"Original","Endpoint phase","Circle trap","Rotated cross trap",
        "Orbit winding","Feedback phase","Finite period","Orbit radius"};
    return i>=0&&i<ColouringCount?names[i]:"Invalid colouring";
}
inline bool finitePoint(C z){return std::isfinite(z.real())&&std::isfinite(z.imag());}
inline C invalidPoint(){const double q=std::numeric_limits<double>::quiet_NaN();return C(q,q);}
// Equivalent to polar(1,arg(z)) for nonzero finite z, including magnitudes
// whose unscaled norm would overflow or underflow. Used only by colouring.
inline C normalizedDirection(C z){
    const double scale=std::max(std::abs(z.real()),std::abs(z.imag()));
    if(scale==0.)return C{};
    const C scaled=z/scale;
    return scaled/std::sqrt(std::norm(scaled));
}
// A chart maps a display point w to its physical point z. No epsilon repairs
// a chart pole. Chart inverses follow the branch nearest the previous point.
inline C chartForward(C w,int chart) {
    if(!finitePoint(w))return invalidPoint();
    C z;
    switch(chart){
    case 0:return w;
    case 1:if(w==C{})return invalidPoint();z=C(1.,0.)/w;break;
    case 2:z=std::exp(w);break;
    case 3:if(w==C{})return invalidPoint();z=.5*(w+C(1.,0.)/w);break;
    default:return invalidPoint();
    }
    return finitePoint(z)?z:invalidPoint();
}
inline C chartInverse(C z,int chart,C reference=C{}) {
    if(!finitePoint(z))return invalidPoint();
    C w;
    switch(chart){
    case 0:return z;
    case 1:if(z==C{})return invalidPoint();w=C(1.,0.)/z;break;
    case 2:{
        if(z==C{}||!finitePoint(reference))return invalidPoint();
        w=std::log(z);
        const double tau=2.*projects::pi;
        w+=C(0.,tau*std::round((reference.imag()-w.imag())/tau));
        break;
    }
    case 3:{
        if(!finitePoint(reference))return invalidPoint();
        using LC=std::complex<long double>;
        const LC physical(z.real(),z.imag()),ref(reference.real(),reference.imag());
        LC root;
        if(std::abs(physical)>1.L){
            const LC inverse=LC(1.L,0.L)/physical;
            root=physical*std::sqrt(LC(1.L,0.L)-inverse*inverse);
        }else root=std::sqrt(physical*physical-LC(1.L,0.L));
        const LC a=physical+root,b=physical-root;
        const LC strong=std::abs(a)>=std::abs(b)?a:b;
        if(strong==LC{})return invalidPoint();
        const LC other=LC(1.L,0.L)/strong; // Root product is exactly one.
        const LC chosen=std::abs(strong-ref)<=std::abs(other-ref)?strong:other;
        w=C(static_cast<double>(chosen.real()),static_cast<double>(chosen.imag()));
        break;
    }
    default:return invalidPoint();
    }
    return finitePoint(w)?w:invalidPoint();
}

inline const char* title(int generation) {
    return generation>=1 && generation<=rhfresh::CandidateCount
        ? rhfresh::registry()[generation-1].name : "RH / invalid candidate";
}
inline const char* formula(int generation) {
    return generation>=1 && generation<=rhfresh::CandidateCount
        ? rhfresh::registry()[generation-1].formula : "Invalid candidate";
}

// Only the mathematical response is replaced. The original native runner's
// orbit loop, smoothing, three palettes, finite-tail display, and bands remain.
class CandidateChannel {
    rhfresh::Model engine_;
    int degree_;
    static rhfresh::Settings settings(int generation,int N,double radius,
            double scale,double shift,C strength,double level) {
        auto s=rhfresh::defaultSettings(generation-1);
        s.candidate=generation-1;s.degree=N;s.radius=radius;s.scale=scale;
        s.shift=shift;s.lambda=strength;
        if(!std::isnan(level))s.level=level;
        return s;
    }
public:
    CandidateChannel(int generation,int N,double radius,double scale,double shift,
                     C strength,double level)
        :engine_(settings(generation,N,radius,scale,shift,strength,level)),degree_(N){}
    bool valid()const{return engine_.valid();}
    int degree()const{return degree_;}
    C response(C z)const{return -engine_.correction(z);}
};

enum class Status { Escaped, Unresolved, NumericalLimit, InvalidModel };
// Public orbit reports retain their full finite-tail diagnostics. Pixel
// rendering requests only the observations consumed by the chosen colour.
enum class OrbitDiagnostics { Full, Colour, LambdaColour };

inline const char* statusName(Status status) {
    switch (status) {
    case Status::Escaped:
        return "display cutoff reached (|z|>16; not certified escape)";
    case Status::Unresolved:
        return "unresolved within iteration budget";
    case Status::NumericalLimit:
        return "nonfinite evaluation or denominator singularity (not escape)";
    default:
        return "invalid model or nonfinite parameter";
    }
}

struct Result {
    C endpoint{};
    int iterations = 0;
    double smooth = 0.0;
    double structure = 0.0;
    double phase = 0.0;
    Status status = Status::Unresolved;
    // Finite-tail recurrence only: not a proof of an attracting cycle.
    int detectedPeriod = 0;
    int periodRepeats = 0;
    double periodResidual = std::numeric_limits<double>::infinity();
    double circleTrap = std::numeric_limits<double>::infinity();
    double crossTrap = std::numeric_limits<double>::infinity();
    double winding = 0.;
    double feedbackPhase = 0., feedbackMean = 0., feedbackPeak = 0.;
    double meanLogRadius = 0.;
    int orbitSamples = 0, feedbackSamples = 0, windingSamples = 0;
};

class Model {
    CandidateChannel channel_;

    static bool finite(C z) {
        return std::isfinite(z.real()) && std::isfinite(z.imag());
    }

    static double logRadius(C z) {
        const double scale = std::max(std::abs(z.real()), std::abs(z.imag()));
        if (scale == 0.0) return -std::numeric_limits<double>::infinity();
        const double re = z.real() / scale, im = z.imag() / scale;
        return std::log(scale) + 0.5 * std::log(re * re + im * im);
    }

    static C power(C z,int exponent) {
        if(exponent==2)return z*z; // Preserve the classic arithmetic exactly.
        C value=z;
        for(int k=1;k<exponent;++k)value*=z;
        return value;
    }

    static double logOnePlusRadius(C z) {
        const double r=logRadius(z);
        if(r==-std::numeric_limits<double>::infinity())return 0.;
        return r>36.?r:std::log1p(std::exp(r));
    }

    static Color errorColour(Status status) {
        return status == Status::InvalidModel ? Color{{.42, .035, .055}}
                                              : Color{{.83, .12, .72}};
    }

    static Color colour(const Result& result, int palette, double hue, double bands) {
        if (result.status == Status::InvalidModel ||
            result.status == Status::NumericalLimit)
            return errorColour(result.status);
        hue = std::isfinite(hue) ? hue : 0.0;
        bands = std::isfinite(bands) ? std::clamp(bands, 0.0, 100.0) : 1.0;
        if (result.status == Status::Unresolved)
            return projects::hsv(hue + .61, .48, .019 + .017 * result.structure);
        const double stripes = .5 + .5 * std::cos(
            bands * (.57 * result.smooth + 2.4 * result.structure));
        const double colorHue = palette == 1 ? hue + .66 - .023 * result.smooth
            : palette == 2 ? hue + .65 * result.structure + .18 * result.phase
            : hue + .57 + .30 * result.phase + .012 * result.smooth;
        return projects::hsv(colorHue, .74, .26 + .72 * stripes);
    }

    struct TailPoint { double real,imag; };
    static void detectTailPeriod(Result& result, const std::array<TailPoint, 80>& tail,
                                 bool alternating=false) {
        if (result.status != Status::Unresolved || result.iterations < 32) return;
        const int available = std::min(result.iterations + 1, 80);
        const int last = result.iterations % 80;
        constexpr double tolerance = 1e-7;
        // Compare four entire consecutive blocks, requiring five blocks of
        // stored states.  The smallest passing period is reported.  Nothing
        // is inferred about later iterates, stability, or exact periodicity.
        for (int period = 1; period <= 16 && 5 * period <= available; ++period) {
            if(alternating&&(period&1))continue;
            double residual = 0.0;
            bool repeated = true;
            for (int lag = 0; lag < 4 * period; ++lag) {
                const auto& a = tail[(last + 80 - lag) % 80];
                const auto& b = tail[(last + 80 - lag - period) % 80];
                const C recent(a.real,a.imag), previous(b.real,b.imag);
                const double error = std::abs(recent - previous) /
                    (1.0 + std::max(std::abs(recent), std::abs(previous)));
                if (!std::isfinite(error) || error > tolerance) {
                    repeated = false;
                    break;
                }
                residual = std::max(residual, error);
            }
            if (repeated) {
                result.detectedPeriod = period;
                result.periodRepeats = 4;
                result.periodResidual = residual;
                return;
            }
        }
    }

    static Color lambdaDiagnosticColour(const Result& result, int palette,
                                        double hue, double bands) {
        if (result.status != Status::Unresolved || result.detectedPeriod == 0)
            return colour(result, palette, hue, bands);
        hue = std::isfinite(hue) ? hue : 0.0;
        bands = std::isfinite(bands) ? std::clamp(bands, 0.0, 100.0) : 1.0;
        const double precision = std::clamp(
            (-std::log10(std::max(1e-16, result.periodResidual)) - 7.0) / 7.0,
            0.0, 1.0);
        const double periodHue = palette == 1
            ? .11 + .057 * result.detectedPeriod
            : .52 + .6180339887498948 * result.detectedPeriod;
        const double phaseHue = palette == 2 ? .09 * result.phase : 0.0;
        const double modulation = .94 + .06 * std::cos(
            bands * (2.0 * projects::pi * result.phase + result.structure));
        const double brightness = (.42 + .30 * precision + .08 * result.structure)
                                  * modulation;
        return projects::hsv(hue + periodHue + phaseHue, .65, brightness);
    }

    Color experimentColour(const Result& r,int palette,double hue,double bands,
                           bool lambdaDiagnostic=false)const {
        if(r.status==Status::InvalidModel||r.status==Status::NumericalLimit)
            return errorColour(r.status);
        if(experiment.coloring==0){
            // Default calls the original function with its original inputs.
            // No multiplication by one changes its floating-point path.
            const double adjustedBands=experiment.frequency==1.?bands:bands*experiment.frequency;
            Color value=lambdaDiagnostic?lambdaDiagnosticColour(r,palette,hue,adjustedBands)
                                        :colour(r,palette,hue,adjustedBands);
            if(experiment.exposure!=1.)
                for(double& component:value)component=std::clamp(component*experiment.exposure,0.,1.);
            return value;
        }
        hue=std::isfinite(hue)?hue:0.;
        bands=std::isfinite(bands)?std::clamp(bands,0.,100.):1.;
        const double freq=experiment.frequency;
        const double tau=2.*projects::pi;
        const double offset=palette==1?.07:palette==2?.36:.59;
        const auto finiteValue=[](double x,double fallback=0.){return std::isfinite(x)?x:fallback;};
        double metric=0.,phase=r.phase,saturation=.76,brightness=.7;
        switch(experiment.coloring){
        case 1:
            metric=r.phase;
            phase=r.phase;
            brightness=r.status==Status::Unresolved?.72:.82;
            break;
        case 2:
            metric=-std::log10(std::max(1e-14,finiteValue(r.circleTrap,1e6)));
            phase=.14*metric;
            brightness=.32+.62*std::exp(-freq*finiteValue(r.circleTrap,1e6));
            break;
        case 3:
            metric=-std::log10(std::max(1e-14,finiteValue(r.crossTrap,1e6)));
            phase=.16*metric;
            brightness=.36+.60*std::exp(-freq*finiteValue(r.crossTrap,1e6));
            break;
        case 4:
            metric=r.winding/tau;
            phase=.17*metric+.12*r.phase;
            brightness=.66;
            break;
        case 5:
            metric=r.feedbackMean;
            phase=r.feedbackPhase+.13*metric;
            brightness=.38+.58*(1.-std::exp(-metric));
            break;
        case 6:
            if(r.detectedPeriod>0){
                metric=static_cast<double>(r.detectedPeriod);
                phase=.6180339887498948*metric;
                brightness=.60+.35*std::clamp(
                    (-std::log10(std::max(1e-16,r.periodResidual))-7.)/7.,0.,1.);
            }else{
                metric=r.smooth;
                phase=.017*r.smooth+.1*r.phase;
                brightness=r.status==Status::Unresolved?.15:.62;
                saturation=.48;
            }
            break;
        case 7:
            metric=r.meanLogRadius;
            phase=.19*metric+.1*r.phase;
            brightness=.32+.63*(1.-std::exp(-metric));
            break;
        default:return errorColour(Status::InvalidModel);
        }
        const double stripes=.5+.5*std::cos(bands*freq*(.71*metric+tau*r.phase));
        const double value=std::clamp(experiment.exposure*brightness*(.66+.34*stripes),0.,1.);
        return projects::hsv(hue+offset+freq*phase,saturation,value);
    }

public:
    int generation = 1;
    int N = 16;
    double x = 0.0;
    double b = 1.0;
    double h = 0.0;
    C lambda{.32, 0.0};
    bool validModel = false;
    int degree = 0;
    std::string description;
    std::string statusText;

    double level = 0.;
    Experiment experiment;
    Model(int g, int truncation, double radius, double scale,
          double shift, C strength,
          double targetLevel=std::numeric_limits<double>::quiet_NaN(),
          const Experiment& exp=Experiment{})
        :channel_(g,truncation,radius,scale,shift,strength,targetLevel),
         generation(g),N(truncation),x(radius),b(scale),h(shift),lambda(strength),
         degree(channel_.degree()),description(title(g)),
         level(std::isnan(targetLevel)?rhfresh::defaultSettings(g-1).level:targetLevel),experiment(exp) {
        validModel=g>=1 && g<=rhfresh::CandidateCount && channel_.valid()
            && b>0. && finite(lambda)&&validExperiment(experiment);
        statusText=validModel
            ? "Finite Theta Taylor model; selected map, chart and orbit colouring experiment."
            : "Require a registered candidate, degree4..32, finite radius/scale>0, finite parameters and a valid experiment.";
    }

    Result orbit(C initial, C parameter, int budget,
                 std::vector<C>* trail = nullptr,
                 OrbitDiagnostics diagnostics=OrbitDiagnostics::Full) const {
        return orbitAtLambda(initial, parameter, lambda, budget, trail,diagnostics);
    }

    Result orbitAtLambda(C initial, C parameter, C strength, int budget,
                         std::vector<C>* trail = nullptr,
                         OrbitDiagnostics diagnostics=OrbitDiagnostics::Full) const {
        Result result;
        result.endpoint = initial;
        if (trail) trail->clear();
        if (!validModel || !validExperiment(experiment) || !finite(parameter) || !finite(strength)) {
            result.status = Status::InvalidModel;
            return result;
        }
        if (!finite(initial)) {
            result.status = Status::NumericalLimit;
            return result;
        }

        budget = std::clamp(budget, 0, 1024);
        C z = initial;
        C previous = initial;
        C feedbackDirection{};
        const C trapRotation=experiment.coloring==3?std::polar(1.,-experiment.trapAngle):C(1.,0.);
        const C feedbackRotation=(experiment.map==9||experiment.map==11)&&experiment.twist!=0.
            ?std::polar(1.,experiment.twist):C(1.,0.);
        const bool fullDiagnostics=diagnostics==OrbitDiagnostics::Full;
        const bool needsTail=fullDiagnostics||experiment.coloring==6
            ||(diagnostics==OrbitDiagnostics::LambdaColour&&experiment.coloring==0);
        const bool needsStructure=fullDiagnostics||experiment.coloring==0;
        // Trivial scalar storage is intentionally uninitialized. Period
        // detection reads only completed states, and most pixel colours
        // never initialize or write these 80 entries at all.
        std::array<TailPoint, 80> tail;
        double minimumTrap = std::numeric_limits<double>::infinity();
        for (int n = 0; n <= budget; ++n) {
            result.endpoint = z;
            result.iterations = n;
            result.smooth = static_cast<double>(n);
            if(needsTail)tail[n % 80] = {z.real(),z.imag()};
            const double largestComponent = std::max(std::abs(z.real()), std::abs(z.imag()));
            if (trail && largestComponent < 1e12)
                trail->push_back(z);
            if(needsStructure)minimumTrap = std::min(minimumTrap, std::abs(z - C(-.21, .37)));
            // Only the selected diagnostic accumulates its extra metric.
            // These observe completed states and never evaluate the driver.
            if(experiment.coloring==2){
                const double distance=std::abs(std::abs(z)-experiment.trapRadius);
                if(std::isfinite(distance))result.circleTrap=std::min(result.circleTrap,distance);
            }else if(experiment.coloring==3){
                const C rotated=z*trapRotation;
                const double distance=std::min(std::abs(rotated.real()),std::abs(rotated.imag()));
                if(std::isfinite(distance))result.crossTrap=std::min(result.crossTrap,distance);
            }else if(experiment.coloring==4){
                if(n>0&&z!=C{}&&previous!=C{}){
                    result.winding+=std::remainder(std::arg(z)-std::arg(previous),2.*projects::pi);
                    ++result.windingSamples;
                }
            }else if(experiment.coloring==7){
                result.meanLogRadius+=logOnePlusRadius(z);
                ++result.orbitSamples;
            }

            // Inputs and every completed update are finite here.  The norm
            // is evaluated only when its components are small.  The smooth
            // colouring log remains representable for any finite endpoint.
            if (largestComponent > 16.0 || std::norm(z) > 256.0) {
                const double logR = logRadius(z);
                result.status = Status::Escaped;
                result.smooth = std::max(0.0, n + 1.0 -
                    std::log(logR / std::log(16.0)) / std::log(2.0));
                break;
            }
            if (n == budget) break;

            C correction{};
            // Preserve the exact unperturbed quadratic map even at a pole
            // of the unused channel: never compute zero times nonfinite.
            C carrier=z;
            if(experiment.map==2||experiment.map==6||(experiment.map==7&&(n&1)))
                carrier=std::conj(z);
            else if(experiment.map==10||experiment.map==11){
                // Exact endpoint branches retain the original holomorphic
                // and conjugate carriers at the blend limits.
                if(experiment.blend==1.)carrier=std::conj(z);
                else if(experiment.blend!=0.)carrier=C(z.real(),(1.-2.*experiment.blend)*z.imag());
            }
            if (strength != C{}&&experiment.gain!=0.) {
                correction = channel_.response(carrier);
                if (!finite(correction)) {
                    result.status = Status::NumericalLimit;
                    break;
                }
            }
            // Keep the original real-strength arithmetic when Im(lambda)=0.
            C feedback = strength.imag() == 0.0
                ? strength.real() * correction : strength * correction;
            if(experiment.gain!=1.)feedback*=experiment.gain;
            if(experiment.map==3&&(n&1))feedback=-feedback;
            if((experiment.map==9||experiment.map==11)&&experiment.twist!=0.)
                feedback*=feedbackRotation;
            if(experiment.coloring==5&&finite(feedback)){
                const double feedbackLog=logOnePlusRadius(feedback);
                result.feedbackMean+=feedbackLog;
                result.feedbackPeak=std::max(result.feedbackPeak,feedbackLog);
                if(feedback!=C{})feedbackDirection+=normalizedDirection(feedback);
                ++result.feedbackSamples;
            }
            C Q;
            if(experiment.map==1){
                Q=z+parameter+feedback;
            }else{
                const C host=power(carrier,experiment.power)+parameter;
                if((experiment.map>=5&&experiment.map<=9)||experiment.map==11){
                    if(feedback==C{})Q=host;
                    else{
                        const C denominator=C(1.,0.)-feedback;
                        if(denominator==C{}){
                            result.status=Status::NumericalLimit;
                            break;
                        }
                        Q=host*((C(1.,0.)+feedback)/denominator);
                    }
                    if(experiment.map==8&&experiment.memory!=0.)Q+=experiment.memory*(z-previous);
                }else{
                    Q=host+feedback;
                    if(experiment.map==4)Q+=experiment.memory*(z-previous);
                }
            }
            // The branch avoids arithmetic by one and subtraction/addition
            // of z in the original, unrelaxed classic path.
            const C next=experiment.relaxation==1.?Q:z+experiment.relaxation*(Q-z);
            if (!finite(next)) {
                result.status = Status::NumericalLimit;
                break;
            }
            previous=z;
            z = next;
        }
        if(needsStructure)result.structure = std::exp(-2.5 * minimumTrap);
        result.phase = .5 + std::arg(result.endpoint) / (2.0 * projects::pi);
        result.phase -= std::floor(result.phase);
        if(result.feedbackSamples>0){
            result.feedbackMean/=result.feedbackSamples;
            result.feedbackPhase=.5+std::arg(feedbackDirection)/(2.*projects::pi);
            result.feedbackPhase-=std::floor(result.feedbackPhase);
        }
        if(result.orbitSamples>0)result.meanLogRadius/=result.orbitSamples;
        if(needsTail)detectTailPeriod(result, tail,experiment.map==3||experiment.map==7);
        return result;
    }

    Color juliaColour(C initial, C parameter, int budget, int palette,
                      double hue, double bands) const {
        return experimentColour(orbit(chartForward(initial,experiment.chart), parameter, budget,
                                      nullptr,OrbitDiagnostics::Colour),
                                palette, hue, bands);
    }

    Color atlasColour(C parameter, int budget, int palette,
                      double hue, double bands) const {
        return experimentColour(orbit(C{}, parameter, budget,nullptr,OrbitDiagnostics::Colour),
                                palette, hue, bands);
    }

    Color lambdaColour(C pixelLambda, C fixedParameter, int budget, int palette,
                       double hue, double bands) const {
        return experimentColour(
            orbitAtLambda(C{}, fixedParameter, pixelLambda, budget,nullptr,OrbitDiagnostics::LambdaColour),
            palette, hue, bands,true);
    }
};

inline std::shared_ptr<const Model> model(int generation,int N=16,
        double radius=.6,double b=.9,double h=.28,C lambda=C(.45,0.),
        double level=std::numeric_limits<double>::quiet_NaN(),
        const Experiment& exp=Experiment{}) {
    static std::array<std::shared_ptr<const Model>,rhfresh::CandidateCount> cache{};
    if(generation<1 || generation>rhfresh::CandidateCount)
        return std::make_shared<const Model>(generation,N,radius,b,h,lambda,level,exp);
    const double target=std::isnan(level)?rhfresh::defaultSettings(generation-1).level:level;
    auto& slot=cache[static_cast<std::size_t>(generation-1)];
    auto cached=std::atomic_load_explicit(&slot,std::memory_order_acquire);
    if(cached&&cached->N==N&&cached->x==radius&&cached->b==b&&cached->h==h
       &&cached->lambda==lambda&&cached->level==target&&sameExperiment(cached->experiment,exp))return cached;
    auto value=std::make_shared<const Model>(generation,N,radius,b,h,lambda,target,exp);
    std::atomic_store_explicit(&slot,value,std::memory_order_release);
    return value;
}

} // namespace rhmaster
