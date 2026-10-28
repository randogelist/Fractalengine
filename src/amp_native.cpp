#ifndef UNICODE
#define UNICODE
#endif
#ifndef _UNICODE
#define _UNICODE
#endif
#include "amp_gpu.hpp"
#include "catalogue.hpp"
#include "self_tests.hpp"
#include "fresh_classic_tests.hpp"
#include "lab_tests.hpp"
#include "gpu_tests.hpp"
#include "focus_presets.hpp"
#include "focus_tests.hpp"
#include "cpu_optimisation_tests.hpp"
#include "render_worker_tests.hpp"
#include "performance_gpu_tests.hpp"
#include "performance_benchmark.hpp"
#include "engine_probe.hpp"
#include "engine_controls_tests.hpp"
#include "stable_curvature_tests.hpp"
#include "stable_curvature_gpu_tests.hpp"
#include "background_gpu.hpp"
#include "rh_branch_extension_tests.hpp"
#include "rh_branch_gpu_tests.hpp"
#include "new_structure_tests.hpp"
#include "new_structure_gpu_tests.hpp"
#include "frk_extra_tests.hpp"
#include "qqrh_hilbert_tests.hpp"
#include "qqrh_hilbert_gpu_tests.hpp"
#include <chrono>
#include <commdlg.h>
#include <condition_variable>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <shellapi.h>
#include <windowsx.h>

using namespace amp;
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
static constexpr UINT FRAME = WM_APP + 1, NOTICE = WM_APP + 2, RENDER_DONE = WM_APP + 3,
    GPU_READY = WM_APP + 4;
static constexpr UINT_PTR INTERACTION_TIMER = 0x5248;
static const wchar_t *title = L"FRAKTALFAKTUR / RH Lab GPU 1.8.6";
static std::string num(double x, int precision = 4) {
    std::ostringstream o;
    o << std::setprecision(precision) << std::fixed << x;
    return o.str();
}
static std::string complexText(C z) {
    return num(z.real()) + (z.imag() < 0 ? " - " : " + ") + num(std::abs(z.imag())) + "i";
}
static const char *masterTitle(int generation) {
    return rhfresh::registry()[std::clamp(generation,1,rhfresh::CandidateCount)-1].name;
}
static void masterModelMetadata(std::ostream &f,const rhmaster::Model *model,C c) {
    if(!model) {f << "RH fresh candidate model: unavailable\n";return;}
    f << "Model: " << masterTitle(model->generation)
      << "\nModel status: " << model->statusText
      << "\nTheta Taylor truncation N: " << model->N
      << "\nNormalized radius R: " << model->x
      << "\nPhysical angular radius r: " << 200*model->x
      << "\nLevel: " << model->level
      << "\nRadius scale b: " << model->b
      << "\nShift h: " << model->h
      << "\nFeedback lambda: " << model->lambda
      << "\nFeedback lambda.real: " << model->lambda.real()
      << "\nFeedback lambda.imag: " << model->lambda.imag()
      << "\nTheta carrier degree: " << model->degree
      << "\nSelected c: " << c
      << "\nCarrier and feedback: " << rhmaster::formula(model->generation)
      << "\nCoordinate: y=b*z+h; T(y)=Theta_N(200y)/Theta(0)."
      << "\nClassic reference update: z <- z^2+c-lambda*F/(b*dF/dy)"
      << "\nMagnitude cutoff |z|>16 is a display cutoff, not a proved escape radius."
      << "\nThe map uses a finite Theta Taylor approximation and rational feedback.\n"
      << "The c portrait follows the single fixed seed z0=0 for each c at the selected lambda.\n"
      << "The lambda portrait follows the single fixed seed z0=0 for each complex lambda at the selected c.\n"
      << "Lambda=1 is a coordinate reference, not an assumed transition.\n"
      << "Detected finite periods are numerical recurrence observations, not certified cycles.\n"
      << "Finite iteration limits and numerical limits are separate unresolved outcomes.\n"
      << "Neither Julia membership, connectedness nor an RH certificate is inferred.\n";
    const auto &e=model->experiment;
    f << "Lab map: " << e.map << " / " << rhmaster::mapName(e.map)
      << "\nLab chart: " << e.chart << " / " << rhmaster::chartName(e.chart)
      << "\nLab coloring: " << e.coloring << " / " << rhmaster::colouringName(e.coloring)
      << "\nLab host power: " << e.power
      << "\nLab feedback gain: " << e.gain
      << "\nLab relaxation: " << e.relaxation
      << "\nLab memory: " << e.memory
      << "\nLab twist (radians; maps 9 / 11): " << e.twist
      << "\nLab conjugate blend (maps 10 / 11): " << e.blend
      << "\nLab trap radius: " << e.trapRadius
      << "\nLab trap angle: " << e.trapAngle
      << "\nLab exposure: " << e.exposure
      << "\nLab color frequency: " << e.frequency
      << "\nCharts transform Julia pixels only; c and lambda parameter planes retain their original coordinates.\n";
}
static std::string masterOrbitSummary(const rhmaster::Result &result) {
    std::string summary=rhmaster::statusName(result.status);
    if(result.detectedPeriod>0)
        summary+=" / finite-tail near-period "+std::to_string(result.detectedPeriod)+
            " ("+std::to_string(result.periodRepeats)+" repeats; residual="+num(result.periodResidual,10)+")";
    return summary;
}
static void directQqrhModelMetadata(std::ostream &f,const qqrhdirect::Model *model,C c) {
    if(!model) {f << "Direct QQRH model: unavailable\n";return;}
    f << "Model: " << qqrhTitle(model->generation) << " / direct"
      << "\nModel status: " << model->statusText
      << "\nModel description: " << model->description
      << "\nCarrier polynomial degree (not full rational-map degree): " << model->degree
      << "\nActual carrier and direct update: " << qqrhDirectFormula(model->generation)
      << "\nTheta Taylor truncation N: " << model->N
      << "\nAngular x: " << model->x << (model->generation==1 || model->generation==6 ? " (unused by curvature)" : "")
      << "\nRadius scale b: " << model->b << "\nShift h: " << model->h
      << "\nFixed complex eta: " << model->eta
      << "\nFixed eta.real: " << model->eta.real() << "\nFixed eta.imag: " << model->eta.imag()
      << "\nMoving-boundary slope m: " << model->slope << (model->generation==5 ? "" : " (unused by this channel)")
      << "\nSelected c: " << c
      << "\nThe c portrait follows z0=0 at each pixel c; the Julia view follows pixel z0 at selected c."
      << "\nMagnitude |z|>16 is a display cutoff, not a proved escape radius."
      << "\nDetected finite periods are numerical recurrence observations, not certified stable cycles."
      << "\nFinite iteration limits and numerical limits are distinct unresolved outcomes."
      << "\nReal-variable positivity certificates do not certify these complex orbits, Julia membership or RH.\n";
    if(model->generation==6)
        f << "Scale kappa: " << model->scale
          << "\nScale-ratio curvature: K_k(y)=k*C_H(k*y)-C_H(y), C_H(y)=-d_y[y*H'(y)/H(y)], H(y)=T(-y)."
          << "\nDirect update: Q=c+eta*[K_k(b*z+h)-K_k(h)]."
          << "\nKappa=1 gives exactly Q=c; kappa=2 is the dyadic case.\n";
}
static const char *legacyMetricTitle(int metric) {
    if(metric==5)return "lyap / finite separation";
    if(metric==7)return "twin / finite separation";
    return bkqrlambda::metricName(metric);
}
static void legacyLambdaMetadata(std::ostream &f,const State &s) {
    const auto settings=bkqrSettings(s);
    if(s.family==23) {
        f << "Model: BKQR lambda wrapper connected to " << qqrhTitle(s.generation)
          << "\nActual carrier and host formula: " << qqrhFormula(s.generation)
          << "\nInner host strength eta: " << s.alpha
          << "\nInner eta.real: " << s.alpha.real() << "\nInner eta.imag: " << s.alpha.imag()
          << "\nOuter lambda: pixel coordinate, independent of eta"
          << "\nTheta Taylor truncation N: " << (8<<s.method)
          << "\nAngular x: " << s.seed.real() << (s.generation==1 ? " (unused by curvature host)" : "")
          << "\nRadius scale b: " << s.seed.imag() << "\nShift h: " << s.spectralD
          << "\nMoving-boundary slope m: " << s.qqrhSlope << (s.generation==5 ? "" : " (unused by this channel)")
          << "\nHost convention: " << (s.generation<=2 ? "Q_j=c+eta*(response_j(b*z+h)-response_j(h)); eta=0 gives the constant host c."
                                                       : "Q_j=z^2+c+eta*Phi_j(z); eta=0 gives the quadratic host z^2+c.")
          << "\nAt outer lambda=0 the update is Q_j. Forced SUB at lambda=1 is identity; GREEDY freezes only on its selected SUB branch for finite Q_j."
          << "\nReal-variable positivity certificates do not certify the complex host orbits or the sampled lambda portrait."
          << "\nThese truncated numerical hosts and ensemble diagnostics do not establish RH.\n";
    } else if(s.family==22) {
        f << "Model: BKQR lambda wrapper connected to " << masterTitle(s.generation)
          << "\nHost Q_j(z): z^2+c+eta*response_j(z)"
          << "\nInner host strength eta: " << s.alpha
          << "\nInner eta.real: " << s.alpha.real() << "\nInner eta.imag: " << s.alpha.imag()
          << "\nOuter lambda: pixel coordinate, independent of eta"
          << "\nTheta Taylor truncation N: " << (8<<s.method)
          << "\nAngular x: " << (s.generation==3 ? -1. : s.seed.real())
          << "\nRadius scale b: " << s.seed.imag() << "\nShift h: " << s.spectralD
          << "\nResponse: -F/(b*dF/dy), y=b*z+h; removable factors cancelled polynomially."
          << "\nBoundary carrier: B_N(y)=A_N(y,-1)."
          << "\nCarrier F: " << (s.generation==1 ? "A_N(y,x)/y" : s.generation==2 ? "(A_N(y,x)-B_N(y))/(y^2*(x+1))" : "B_N(y)/y, B_N(y)=A_N(y,-1)")
          << "\nAt outer lambda=0 the update is Q_j. Forced SUB at lambda=1 is the identity (including its continuous pole extension)."
          << "\nGREEDY at lambda=1 freezes only where Q_j is finite and the SUB branch is chosen; it is not universal collapse."
          << "\nThese finite RH-host approximations and wrapper diagnostics are not RH certificates.\n";
    } else {
        f << "Model: original BKQR lambda map\n"
          << "Host: " << bkqrlambda::hostName(settings.host) << "\nHost power: " << settings.power << '\n';
    }
    f << "Branch policy: " << bkqrlambda::branchName(settings.branch)
      << "\nLambda coordinates: " << bkqrlambda::couplingName(settings.coupling)
      << "\nHost parameter c: " << settings.c
      << "\nStep budget: " << settings.iterations
      << "\nSampling quality: " << settings.quality
      << "\nPrecision preset: " << bkqrlambda::precisionName(settings.precision)
      << " (grid/seed budget; CPU double arithmetic)"
      << "\nRequested grid side: " << bkqrlambda::gridSize(settings)
      << "\nInitial seeds per grid point: " << bkqrlambda::seedCount(settings)
      << "\nMetric: " << bkqrlambda::metricName(settings.metric)
      << "\nUpdate: z_next=z+(1+sigma*H(lambda))*(Q(z)-z)."
      << "\nCoordinate H(lambda): " << (settings.coupling==0 ? "lambda" : settings.coupling==1 ? "4*lambda/(1+lambda)^2" : "2*lambda/(1+lambda)")
      << "\nBranch sigma: " << (settings.branch==0 ? "-1 if Re(Q(z)*conj(Q(z)-z))>0; otherwise +1" : settings.branch==1 ? "-1" : "+1")
      << "\nCoordinate display: centered square; equal x/y scales; lambda plane only."
      << "\nCandidate is a gradient of sampled statistics, not a certified phase boundary."
      << "\nLyap and twin use the same finite, unrenormalized shadow-separation metric, not a rigorous Lyapunov exponent."
      << "\nEscape retains the legacy finite/numerical stop test; survival is the normalized executed-step count."
      << "\nFinite ensemble statistics, step budgets, and numerical limits do not certify infinite-time behavior.\n";
}
static const char *riemannTitle(int generation) {
    static const char *names[]={"", "RH / finite moment witness nodes", "RH / cubic pullback basins",
        "RH / reciprocal cubic sheets", "RH / orthogonal pencil", "RH / mirrored moment kernel",
        "RH / Jacobi closure", "RH / Joukowski pencil"};
    return names[std::clamp(generation,1,7)];
}
static void riemannModelMetadata(std::ostream &f,const riemann::Model *model) {
    if(!model) {f << "RH model: unavailable\n";return;}
    const auto previous=f.precision();
    f << std::setprecision(std::numeric_limits<long double>::max_digits10)
      << "Model: " << model->description << "\nModel status: " << model->statusText
      << "\nCoefficient provenance: " << model->coefficientOrigin
      << "\nPolynomial degree: " << model->degree << "\nMaximum relative moment residual: " << model->maxMomentResidual << '\n';
    auto values=[&](const char *name,const auto &list) {
        f << name << ':';
        for(const auto &v:list) f << ' ' << v;
        f << '\n';
    };
    values("Frozen source moments",model->moments);
    values("Witness nodes q",model->nodes);
    values("Moment weights",model->weights);
    values("Ordinary weights",model->ordinaryWeights);
    if(model->generation>=4) {
        f << "Standardization mean mu: " << model->standardizedMean
          << "\nStandardization scale sigma: " << model->standardizedScale
          << "\nStandardized skew b=E[t^3]: " << model->standardizedSkew
          << "\nStandardized H h=E[t^4]-b^2-1: " << model->standardizedH
          << "\nCoordinate: " << (model->generation==7 ? "Joukowski variable z, t=(z+1/z)/2, t=(q-mu)/sigma"
                                                       : "z=(q-mu)/sigma")
          << "\nPolynomial inputs: five fixed moments c0..c4; independent of the witness anchor.\n";
        if(model->generation==6)
            f << "Jacobi closure: terminal diagonal is set to zero; it is an explicit finite closure choice.\n";
        if(model->generation==7)
            f << "The cleared Joukowski expression is a degree-eight polynomial, including at z=0.\n";
    }
    f << "Residual=max_n |sum_j W_j q_j^(n+1)-c_n|/|c_n| for the stored finite moments.\n"
      << "See RIEMANN_WITNESS_MATHEMATICS.txt for source and construction details.\n"
      << "Coefficients are constructed numerically from the bundled finite RH-branch moment data.\n"
      << "The finite witness and its basin geometry do not establish the infinite RH binding.\n";
    f.precision(previous);
}
static void quotientModelMetadata(std::ostream &f,const momentquotient::Model *model,C c) {
    if(!model) {f << "RH quotient model: unavailable\n";return;}
    f << "Model: " << model->description << "\nModel status: " << model->statusText
      << "\nIteration: f_c(z)=lambda*(z-b-1/z)+c"
      << "\nb=E[t^3]: " << model->b << "\nh=E[t^4]-b^2-1: " << model->h
      << "\nlambda=1/sqrt(h): " << model->lambda
      << "\nInfinity multiplier (local coordinate 1/z): " << model->infinityMultiplier
      << "\nCritical points: +i, -i\nPole: 0 -> infinity\nInfinity: fixed"
      << "\nSelected c: " << c << "\nEscape radius for selected c: " << model->escapeRadius(c) << '\n';
    f << "BEGIN FIXED MOMENT SOURCE (not the iterated map)\n";
    riemannModelMetadata(f,model->source.get());
    f << "END FIXED MOMENT SOURCE\n";
}
static void tripleModelMetadata(std::ostream &f,const momenttriple::Model *model,C c) {
    if(!model) {f << "RH triple model: unavailable\n";return;}
    f << "Model: " << model->description << "\nModel status: " << model->statusText
      << "\nIteration: " << momenttriple::formula(model->generation)
      << "\nb=E[t^3]: " << model->b << "\nh=E[t^4]-b^2-1: " << model->h
      << "\nA=1+h+b^2/3: " << model->A << "\na=1+b^2/4: " << model->a
      << "\neta=sqrt(h): " << model->eta
      << "\nDegree (absolute for anti-map): " << model->degree
      << "\nPower at infinity: " << model->infinityPower
      << "\nSelected c: " << c << "\nEscape radius for selected c: " << model->escapeRadius(c)
      << "\nStep unit: " << (model->generation==2 ? "one full return / two alternating quadratic updates" : "one map update")
      << "\nFinite critical points, ordered for atlas bits:\n";
    const auto starts=model->criticalPoints(c);
    for(size_t k=0;k<starts.size();++k)f << k << ": " << starts[k] << '\n';
    f << "Infinity: fixed and superattracting; always in its own basin.\n";
    if(model->generation==3)f << "Pole: 0 -> infinity.\n";
    if(model->generation==1)f << "Antiholomorphic update; Julia set understood via its holomorphic second iterate.\n";
    f << "BEGIN FIXED MOMENT SOURCE (not the iterated map)\n";
    riemannModelMetadata(f,model->source.get());
    f << "END FIXED MOMENT SOURCE\n";
}
static const char *projectTitle(int generation) {
    static const char *names[]={"", "GRH / completed-pair polynomial model", "PI / actual-digit discrepancy basins",
        "NSE / diffusing spectral-energy basins", "COLLATZ / entire two-step dynamics", "P_c(d,q) / critical-dimension basins"};
    return names[std::clamp(generation,1,5)];
}
static std::string periodText(int p) { return p ? std::to_string(p) : "unresolved"; }
static std::string stamp() {
    SYSTEMTIME t;
    GetLocalTime(&t);
    char s[64];
    std::snprintf(s, sizeof s, "%04u%02u%02u_%02u%02u%02u_%03u", t.wYear, t.wMonth, t.wDay, t.wHour,
                  t.wMinute, t.wSecond, t.wMilliseconds);
    return s;
}
static void metadata(const fs::path &p, const State &s, int w, int h, int samples,
                     const std::string &backend, const quantum::Stats *stats = nullptr,
                     const spherical::Stats *sphereStats = nullptr,
                     const spin::Stats *spinStats = nullptr,
                     const projects::Model *projectModel = nullptr,
                     const riemann::Model *riemannModel = nullptr,
                     const momentcubic::Model *cubicModel = nullptr,
                     const momentquotient::Model *quotientModel = nullptr,
                     const momenttriple::Model *tripleModel = nullptr,
                     const rhmaster::Model *masterModel = nullptr,
                     const qqrhdirect::Model *qqrhDirectModel = nullptr) {
    if(s.family==25) {
        std::ofstream f(fs::path(p.wstring()+L".txt"));
        if(!f)throw std::runtime_error("Cannot write new fractal metadata");
        const auto &camera=s.cameras[s.plane];
        const auto &entry=frk_extra::extra_catalogue[std::clamp(s.generation,1,frk_extra::driver_count)-1];
        f << std::setprecision(17) << "FRAKTALFAKTUR RH Lab GPU 1.8.6 / integrated additions\n"
          << "Driver: " << entry.driver_id << " / " << entry.label << '\n'
          << formula(s) << "\nPlane: " << (s.plane==0 ? "parameter c" : "initial state z")
          << "\nCamera: " << camera.x << ' ' << camera.y << ' ' << camera.span
          << "\nSelected c: " << s.c << "\nReal lambda: " << s.alpha.real()
          << "\nSize: " << w << 'x' << h << "\nSamples per axis: " << samples
          << "\nIterations: " << s.iterations << "\nPalette: " << s.palette
          << "\nBackend: " << backend
          << "\nNonescape and finite-time root convergence are numerical observations; no theorem is certified.\n";
        if(s.generation>=9) f << "Finite Hardy model: N=12 (13 components)\nScaled radius R: "
          << s.alpha.real() << "\nPhysical Theta radius r=200R: " << 200.*s.alpha.real()
          << "\nData: original theta_coefficients.hpp, floating-point b_n=200^n*a_n/a_0\n"
          << "Model: normalized v_n=b_n*R^n; c_k=<B^k*v,v>. No full-function tail bound.\n";
        if(!f)throw std::runtime_error("New fractal metadata write failed");
        writeView(p.wstring()+L".view",s);return;
    }
    if(s.family==24) {
        std::ofstream f(fs::path(p.wstring()+L".txt"));
        if(!f)throw std::runtime_error("Cannot write direct QQRH metadata");
        const auto &v=s.cameras[s.plane];
        f << std::setprecision(17) << "FRAKTALFAKTUR Native Catalogue 0.24.0\n"
          << qqrhTitle(s.generation) << " / direct\n" << qqrhDirectFormula(s.generation)
          << "\nPlane: " << (s.plane==0 ? "parameter c / fixed seed z0=0 / fixed complex eta" : "initial state z / selected c / fixed complex eta")
          << "\nCamera: " << v.x << ' ' << v.y << ' ' << v.span
          << "\nSize: " << w << 'x' << h << "\nSamples per axis: " << samples
          << "\nIteration budget: " << s.iterations << "\nBackend: " << backend
          << "\nPalette: " << s.palette << "\nHue/bands: " << s.hue << ' ' << s.bands << '\n';
        directQqrhModelMetadata(f,qqrhDirectModel,s.c);
        if(!f)throw std::runtime_error("Direct QQRH metadata write failed");
        writeView(p.wstring()+L".view",s);return;
    }
    if(s.family==21 || s.family==22 || s.family==23) {
        std::ofstream f(fs::path(p.wstring()+L".txt"));
        if(!f) throw std::runtime_error("Cannot write BKQR lambda-map metadata");
        const auto &v=s.cameras[2];
        f << std::setprecision(17) << "FRAKTALFAKTUR Native Catalogue 0.24.0\n"
          << formula(s) << "\nPlane: " << (s.family==23 ? "BKQR outer complex lambda / QQ-RH host" : s.family==22 ? "BKQR outer complex lambda / RH host" : "original BKQR complex lambda map")
          << "\nCamera centre/span: " << v.x << ' ' << v.y << ' ' << v.span
          << "\nLambda real bounds: " << v.x-v.span/2 << ' ' << v.x+v.span/2
          << "\nLambda imaginary bounds: " << v.y-v.span/2 << ' ' << v.y+v.span/2
          << "\nCanvas size: " << w << 'x' << h << "\nSquare display side: " << std::min(w,h)
          << "\nCanvas projection: nearest finite-grid sample; a larger canvas does not add grid samples."
          << "\nBackend: " << backend << '\n';
        legacyLambdaMetadata(f,s);
        if(!f) throw std::runtime_error("BKQR lambda-map metadata write failed");
        writeView(p.wstring()+L".view",s);return;
    }
    if(s.family==20) {
        std::ofstream f(fs::path(p.wstring()+L".txt"));
        if(!f) throw std::runtime_error("Cannot write RH fresh candidate metadata");
        const auto &v=s.cameras[s.plane];
        f << std::setprecision(17) << "FRAKTALFAKTUR RH Lab GPU 1.8.0\n"
          << masterTitle(s.generation) << '\n' << formula(s)
          << "\nPlane: " << (s.plane==0 ? "parameter c / fixed seed z0=0 / selected complex lambda" :
                               s.plane==2 ? "complex lambda / fixed seed z0=0 / selected c" : "initial state z / selected c and complex lambda")
          << "\nCamera: " << v.x << ' ' << v.y << ' ' << v.span
          << "\nSize: " << w << 'x' << h << "\nSamples per axis: " << samples
          << "\nIteration budget: " << s.iterations << "\nBackend: " << backend
          << "\nPalette: " << s.palette << "\nHue/bands: " << s.hue << ' ' << s.bands << '\n';
        masterModelMetadata(f,masterModel,s.c);
        if(!f)throw std::runtime_error("RH fresh candidate metadata write failed");
        writeView(p.wstring()+L".view",s);return;
    }
    if(s.family==19) {
        std::ofstream f(fs::path(p.wstring()+L".txt"));
        if(!f) throw std::runtime_error("Cannot write RH three-way metadata");
        const auto &v=s.cameras[s.plane];
        f << std::setprecision(17) << "FRAKTALFAKTUR Native Catalogue 0.24.0\n"
          << momenttriple::title(s.generation) << '\n' << formula(s)
          << "\nPlane: " << (s.plane==0 ? "parameter c / all finite critical orbits" : "initial state z / infinity-basin diagnostics")
          << "\nCamera: " << v.x << ' ' << v.y << ' ' << v.span
          << "\nSize: " << w << 'x' << h << "\nSamples per axis: " << samples
          << "\nIteration budget: " << s.iterations << "\nBackend: " << backend
          << "\nPalette: " << s.palette << "\nHue/bands: " << s.hue << ' ' << s.bands << '\n';
        tripleModelMetadata(f,tripleModel,s.c);
        f << "Each atlas pixel uses its own c, critical starts and escape radius.\n"
          << "Palette 0 hue encodes the escaped-critical subset; dark means none observed escaping.\n"
          << "Finite nonescape is unresolved; no exact Julia membership or connectedness is certified.\n";
        if(!f)throw std::runtime_error("RH three-way metadata write failed");
        writeView(p.wstring()+L".view",s);return;
    }
    if(s.family==18) {
        std::ofstream f(fs::path(p.wstring()+L".txt"));
        if(!f) throw std::runtime_error("Cannot write RH quotient metadata");
        const auto &v=s.cameras[s.plane];
        f << std::setprecision(17) << "FRAKTALFAKTUR Native Catalogue 0.24.0\n"
          << "RH / Jacobi quotient Julia atlas\n" << formula(s)
          << "\nPlane: " << (s.plane==0 ? "parameter c / critical escape atlas" : "initial state z / infinity basin")
          << "\nCamera: " << v.x << ' ' << v.y << ' ' << v.span
          << "\nSize: " << w << 'x' << h << "\nSamples per axis: " << samples
          << "\nIteration budget: " << s.iterations << "\nBackend: " << backend
          << "\nPalette: " << s.palette << "\nHue/bands: " << s.hue << ' ' << s.bands << '\n';
        quotientModelMetadata(f,quotientModel,s.c);
        f << "The atlas follows both critical points; each pixel has its own c and escape radius.\n"
          << "Homogeneous coordinates include pole 0 and the fixed point at infinity.\n"
          << "The coloured escape region approximates the infinity basin; finite nonescape remains unresolved.\n"
          << "No rational connectedness conclusion follows from this finite classification.\n";
        if(!f) throw std::runtime_error("RH quotient metadata write failed");
        writeView(p.wstring()+L".view",s);
        return;
    }
    if(s.family==17) {
        std::ofstream f(fs::path(p.wstring()+L".txt"));
        if(!f) throw std::runtime_error("Cannot write RH cubic metadata");
        const auto &v=s.cameras[s.plane];
        f << std::setprecision(17) << "FRAKTALFAKTUR Native Catalogue 0.24.0\n"
          << "RH / Jacobi cubic Julia atlas\n" << formula(s)
          << "\nPlane: " << (s.plane==0 ? "parameter c / both critical orbits" : "initial state z / selected c")
          << "\nSelected c: " << s.c << "\nCamera: " << v.x << ' ' << v.y << ' ' << v.span
          << "\nSize: " << w << 'x' << h << "\nSamples per axis: " << samples
          << "\nIteration budget: " << s.iterations << "\nBackend: " << backend
          << "\nPalette: " << s.palette << "\nHue/bands: " << s.hue << ' ' << s.bands << '\n';
        if(cubicModel) {
            f << "Model: " << cubicModel->description << "\nModel status: " << cubicModel->statusText
              << "\nA=1+h+b^2/3: " << cubicModel->A
              << "\nCritical points: +/- " << cubicModel->critical
              << "\nEscape radius for selected c: " << cubicModel->escapeRadius(s.c) << '\n';
            f << "BEGIN FIXED MOMENT SOURCE (not the iterated map)\n";
            riemannModelMetadata(f,cubicModel->source.get());
            f << "END FIXED MOMENT SOURCE\n";
        }
        f << "Iteration: z <- z^3-A*z+c; no Newton or Halley correction.\n"
          << "The atlas uses both critical points; each pixel has its own c and escape radius.\n"
          << "Escape is a finite floating-point radius test. A budget-limited orbit is unresolved.\n"
          << "No connectedness or infinite-orbit conclusion is certified by nonescape.\n";
        if(!f) throw std::runtime_error("RH cubic metadata write failed");
        writeView(p.wstring()+L".view",s);
        return;
    }
    if(s.family==16) {
        std::ofstream f(fs::path(p.wstring()+L".txt"));
        if(!f) throw std::runtime_error("Cannot write RH witness metadata");
        const auto &v=s.cameras[0];
        f << std::setprecision(17) << "FRAKTALFAKTUR Native Catalogue 0.24.0\n"
          << riemannTitle(s.generation) << '\n' << formula(s)
          << "\nGeneration: " << s.generation << "\nMoment order: H" << s.riemann.order
          << (s.generation>=4 ? "\nH2 anchor (witness diagnostics only): " : "\nH2 prescribed anchor: ") << s.riemann.anchor
          << "\nMethod: " << (s.method ? "Halley" : "Newton") << "\nComplex relaxation: " << s.alpha
          << "\nPlane: complex initial state z\nCamera: " << v.x << ' ' << v.y << ' ' << v.span
          << "\nSize: " << w << 'x' << h << "\nSamples per axis: " << samples
          << "\nIteration budget: " << s.iterations << "\nBackend: " << backend
          << "\nPalette: " << s.palette << "\nHue/bands: " << s.hue << ' ' << s.bands << '\n';
        riemannModelMetadata(f,riemannModel);
        f << "Finite numerical root basins; iteration limits and evaluation cutoffs remain unresolved.\n";
        if(!f) throw std::runtime_error("RH witness metadata write failed");
        writeView(p.wstring()+L".view",s);
        return;
    }
    if (s.family == 15) {
        std::ofstream f(fs::path(p.wstring()+L".txt"));
        if (!f) throw std::runtime_error("Cannot write project-fractal metadata");
        const auto &t=s.project;
        const auto &v=s.cameras[0];
        f << std::setprecision(17) << "FRAKTALFAKTUR Native Catalogue 0.24.0\n"
          << projectTitle(s.generation) << '\n' << formula(s)
          << "\nProject index: " << s.generation << "\nPlane: complex initial state z"
          << "\nCamera: " << v.x << ' ' << v.y << ' ' << v.span
          << "\nSize: " << w << 'x' << h << "\nSamples per axis: " << samples
          << "\nIteration budget: " << s.iterations << "\nBackend: " << backend
          << "\nPalette: " << s.palette << "\nHue/bands: " << s.hue << ' ' << s.bands << '\n';
        if(projectModel) f << "Model: " << projectModel->description << "\nModel status: " << projectModel->statusText
            << "\nPolynomial degree: " << projectModel->degree << "\nPi digit certificate: "
            << (s.generation==2 ? (projectModel->piCertified ? "PASS" : "FAIL") : "not applicable") << '\n';
        if (s.generation==1)
            f << "GRH polynomial degree: " << t.grhDegree << "\nGRH coefficient scale: " << t.grhRadius
              << "\nOffset a: " << s.c << "\nMultiplier b: " << s.alpha
              << "\nFinite polynomial approximation of the completed character pair; dynamics are those of the approximation.\n";
        else if (s.generation==2)
            f << "Pi terms: " << t.piTerms << "\nPi pair offset: " << t.piOffset
              << "\nFixed drift anchor N0=1024; n=N0+pair_offset+t. Coefficients are beta_(2n)-beta_(2n+1)-N0*c(n).\n";
        else if (s.generation==3)
            f << "NSE mode count: " << t.nseModes << "\nNSE nu*t: " << t.nseTime
              << "\nNSE energy slope: " << t.nseSlope
              << "\nExact periodic shear sector; the picture is the Newton iteration of its finite spectral-energy polynomial.\n";
        else if (s.generation==4)
            f << "Entire extension F(z)=z+1/4-(z/2+1/4)cos(pi*z); each displayed step is F(F(z)).\n";
        else
            f << "PCD q=2 probability p: " << t.pcdProbability
              << "\nComplex sigma=2*d-1; Newton basins of the GP-4 quintic.\n";
        f << "Finite numerical classification; iteration limits, poles and magnitude cutoffs are not escape proofs.\n"
          << "The visual dynamics do not certify GRH, pi non-normality, NSE regularity, Collatz, or a physical percolation threshold.\n";
        if (!f) throw std::runtime_error("Project-fractal metadata write failed");
        writeView(p.wstring()+L".view",s);
        return;
    }
    if (s.family == 14) {
        std::ofstream f(fs::path(p.wstring() + L".txt"));
        if (!f) throw std::runtime_error("Cannot write spin dynamics metadata");
        const auto &q=s.quantum;
        const auto &t=s.spin;
        const bool atlas=s.generation==3 && s.plane==0;
        f << std::setprecision(17) << "FRAKTALFAKTUR Native Catalogue 0.24.0\n"
          << "Dissipative kicked top / finite numerical diagnostics\n" << formula(s)
          << "\nGeneration: " << s.generation << "\nTwist k: " << t.twist
          << "\nTilt beta: " << t.tilt << "\nDamping tau: " << t.damping
          << "\nView: " << s.plane
          << (s.generation==3 ? " (0=parameter atlas,1=selected sphere,2=selected equal-area atlas)"
                              : " (0=sphere,1=equal-area atlas,2=front/back hemispheres)")
          << "\nAtlas k range: " << t.kMin << ' ' << t.kMax
          << "\nAtlas tau range: " << t.tauMin << ' ' << t.tauMax
          << "\nGrid power: " << t.gridPower << "\nRequested burn-in: " << (uint64_t(1)<<t.burnPower)
          << "\nLyapunov averaging window: " << (uint64_t(1)<<t.windowPower)
          << "\nInitial-condition seed: " << q.randomSeed
          << "\nCamera yaw/pitch/zoom: " << q.yaw << ' ' << q.pitch << ' ' << q.zoom
          << "\nExposure: " << q.exposure << "\nPalette/hue: " << s.palette << ' ' << s.hue
          << "\nGuides: " << s.method << "\nSize: " << w << 'x' << h << "\nBackend: " << backend
          << "\nRequested base samples: " << (uint64_t(1)<<q.samplePower)
          << "\nHQ: " << (atlas ? "denser parameter grid; see realized grid statistics" : "up to 4x retained points, capped at 2^23") << '\n';
        f << "Export sample factor per axis: " << samples << '\n';
        if (spinStats) {
            f << "Actual discarded burn-in: " << spinStats->burnIn
              << "\nActual averaging window: " << spinStats->window
              << "\nMean finite-time Lyapunov estimate: " << spinStats->meanLyapunov
              << "\nMean log area rate: " << spinStats->meanAreaRate
              << "\nFinite-time Lyapunov range: " << spinStats->minLyapunov << ' ' << spinStats->maxLyapunov << '\n';
            if (spinStats->parameterAtlas)
                f << "Grid width/height: " << spinStats->gridWidth << ' ' << spinStats->gridHeight
                  << "\nParameter cells: " << spinStats->cells
                  << "\nInitial states per cell: " << spinStats->startsPerCell
                  << "\nReported means aggregate parameter cells; there is no single cloud mean.\n";
            else
                f << "Retained samples: " << spinStats->samples << "\nIndependent starts: " << spinStats->chains
                  << "\nMean state: " << spinStats->mean.x << ' ' << spinStats->mean.y << ' ' << spinStats->mean.z
                  << "\nMean state norm: " << spinStats->meanNorm
                  << "\nReported growth values are rolling last-window estimates attached to retained points.\n";
        }
        f << (atlas ? "Parameter atlas uses the same one seeded initial state in every cell; coexisting basins can be missed.\n"
                    : "The cloud is a finite occupation sample after burn-in.\n")
          << "The seed changes initial conditions only; every subsequent map step is deterministic.\n"
          << "Positive finite-time growth is not a proof of chaos or fractal dimension.\n";
        if (!f) throw std::runtime_error("Spin metadata write failed");
        writeView(p.wstring()+L".view",s);
        return;
    }
    if (s.family == 13) {
        std::ofstream f(fs::path(p.wstring() + L".txt"));
        if (!f) throw std::runtime_error("Cannot write spherical fractal metadata");
        const auto &q=s.quantum;
        f << std::setprecision(17) << "FRAKTALFAKTUR Native Catalogue 0.24.0\n"
          << "Spherical fractals / finite numerical sample\n" << formula(s)
          << "\nGeneration: " << s.generation
          << "\nSampling: " << (s.generation==3 ? "uniform inverse branches of the degree-six rational map"
                                                      : "uniform reduced reflection words; immediate repeat excluded")
          << "\nMeasure: " << (s.generation==3 ? "empirical inverse-iteration measure" : "empirical group-word visitation measure")
          << "\nSchottky gap (generation 2 only): " << s.sphere.gap
          << "\nLambda log10 magnitude/phase (generation 3 only): " << s.sphere.lambdaLog << ' ' << s.sphere.lambdaPhase
          << "\nConformal display boost: " << s.sphere.boost
          << "\nCamera yaw/pitch/zoom: " << q.yaw << ' ' << q.pitch << ' ' << q.zoom
          << "\nExposure: " << q.exposure << "\nPalette/hue: " << s.palette << ' ' << s.hue
          << "\nView: " << s.plane << " (0=sphere,1=equal-area atlas,2=front/back hemispheres)"
          << "\nGuides: " << s.method << "\nReplay seed: " << q.randomSeed
          << "\nRequested base samples: " << (uint64_t(1)<<q.samplePower)
          << "\nSize: " << w << 'x' << h << "\nBackend: " << backend
          << "\nHQ sampling multiplier: " << samples*samples << " (capped at 2^23)\n";
        if (sphereStats) {
            f << "Retained samples: " << sphereStats->samples
              << "\nIndependent chains: " << sphereStats->chains
              << "\nDiscarded burn-in per chain: " << sphereStats->burnIn
              << "\nUnboosted sample mean: " << sphereStats->mean.x << ' ' << sphereStats->mean.y << ' ' << sphereStats->mean.z
              << "\nUnboosted sample mean norm: " << sphereStats->meanNorm << '\n';
            const int branches=s.generation==3 ? 6 : 4;
            for (int i=0;i<branches;++i) f << "Retained branch " << i+1 << ": " << sphereStats->branchCounts[i] << '\n';
        }
        f << "Conformal boost changes the displayed coordinates, not the sampled dynamics.\n"
          << "These densities are not Born probabilities or uniform surface-area samples.\n"
          << "Finite sampling does not certify stationarity, dimension or parameter-specific topology.\n";
        if (!f) throw std::runtime_error("Spherical metadata write failed");
        writeView(p.wstring()+L".view",s);
        return;
    }
    if (s.family == 12) {
        std::ofstream f(fs::path(p.wstring() + L".txt"));
        if (!f) throw std::runtime_error("Cannot write measurement metadata");
        const auto &q = s.quantum;
        f << std::setprecision(17) << "FRAKTALFAKTUR Native Catalogue 0.24.0\n"
          << "Quantum measurement histories / finite numerical sample\n" << formula(s)
          << "\nGeometry: " << quantum::geometryName(q.geometry)
          << "\nStrength: " << q.strength << "\nKick radians: " << q.kick
          << "\nKick axis azimuth/elevation: " << q.axisAzimuth << ' ' << q.axisElevation
          << "\nCamera yaw/pitch/zoom: " << q.yaw << ' ' << q.pitch << ' ' << q.zoom
          << "\nExposure: " << q.exposure << "\nPalette/hue: " << s.palette << ' ' << s.hue
          << "\nView: " << s.plane << " (0=sphere,1=equal-area atlas,2=mean/history)"
          << "\nGuides: " << s.method << "\nReplay seed: " << q.randomSeed
          << "\nRequested base samples: " << (uint64_t(1) << q.samplePower)
          << "\nSize: " << w << 'x' << h << "\nBackend: " << backend
          << "\nHQ sampling multiplier: " << samples*samples << " (capped at 2^23)\n";
        if (stats) f << "Retained samples: " << stats->samples << "\nIndependent chains: " << stats->chains
          << "\nDiscarded burn-in per chain: " << stats->burnIn
          << "\nRetained-history mean: " << stats->mean.x << ' ' << stats->mean.y << ' ' << stats->mean.z
          << "\nRetained-history mean norm: " << stats->meanNorm << '\n';
        f << "Exact ensemble factor per step: " << quantum::meanFactor(q.strength)
          << "\nThe mean chart uses an ensemble starting north; retained density samples follow burn-in.\n"
          << "Finite sampling does not certify stationarity, dimension or uniqueness of a limiting measure.\n";
        if (!f) throw std::runtime_error("Measurement metadata write failed");
        writeView(p.wstring()+L".view",s);
        return;
    }
    auto v = s.cameras[s.plane];
    std::ofstream f(fs::path(p.wstring() + L".txt"));
    if (!f)
        throw std::runtime_error("Cannot write render metadata");
    f << std::setprecision(17) << "FRAKTALFAKTUR Native Catalogue 0.24.0\n"
      << formula(s) << "\nFamily: " << familyName(s.family) << "\nGeneration: " << s.generation
      << "\nMethod: " << (s.family == 11 ? "Direct spectral jet iteration" : s.family == 9 ? "Three-panel direct diagnostic" : s.family == 10 ? (s.method ? "Halley + phase feedback" : "Newton + phase feedback") : s.family ? (s.method ? "Halley" : "Newton") : "Direct") << "\n";
    if (s.family == 9)
        f << "Plane: shared x=log2(a), y=logit(r)\nSize: " << w << " x " << h;
    else
        f << "Plane: " << s.plane << (s.family == 11 ? " (0=x-start, 1=c, 2=eta)\nSize: " : " (0=z-start, 1=c, 2=alpha)\nSize: ") << w << " x " << h;
    f << "\nSamples per axis: " << samples << "\nIterations: " << s.iterations << "\nBackend: " << backend
      << "\nc: " << s.c << "\nalpha: " << s.alpha << "\nCamera: " << v.x << ' ' << v.y << ' ' << v.span
      << "\n";
    if (s.family == 11)
        f << "Spectral spacing d: " << s.spectralD << "\nOuter atoms: " << .5-s.spectralD << ' ' << .5+s.spectralD
          << "\neta: " << s.alpha << "\nParameter seed: " << s.seed << "\n";
    if (s.family == 11 && s.generation == 2)
        f << "Middle atom: " << .5+.35*s.spectralD << "\nMiddle weight: " << s.spectralWeight
          << "\nJet order: " << s.spectralOrder << "\nWeights: " << (1-s.spectralWeight)/2 << ' '
          << s.spectralWeight << ' ' << (1-s.spectralWeight)/2 << "\n";
    f << "Finite numerical diagnostics, not a proof of period or fractality.\n"
      << (s.family == 11 ? "Direct rational dynamics; parameter plots track the selected seed. Nonescaped pixels remain unclassified.\n"
          : s.family == 9 ? "Literal observable visualization; floating-point evaluation is not a GRH certificate.\n"
          : s.family == 10 ? "Exploratory phase-feedback dynamics; the added feedback is not a GRH equivalence.\n"
          : s.family ? "Numerical root iteration; cutoffs do not prove divergence.\n"
                     : "Homogeneous sphere iteration; no escape-radius test.\n");
    if (!f)
        throw std::runtime_error("Metadata write failed");
    writeView(p.wstring() + L".view", s);
}
struct Frame {
    uint64_t generation = 0;
    int w = 0, h = 0, targetW = 0, targetH = 0, stride = 1;
    double ms = 0, inputMs = -1;
    rhgpu::RenderMetrics gpuMetrics;
    Clock::time_point inputAt{};
    std::vector<uint32_t> pixels;
    std::string backend;
    State state;
    quantum::Stats quantumStats;
    spherical::Stats sphereStats;
    spin::Stats spinStats;
    std::shared_ptr<const projects::Model> projectModel;
    std::shared_ptr<const riemann::Model> riemannModel;
    std::shared_ptr<const momentcubic::Model> cubicModel;
    std::shared_ptr<const momentquotient::Model> quotientModel;
    std::shared_ptr<const momenttriple::Model> tripleModel;
    std::shared_ptr<const rhmaster::Model> masterModel;
    std::shared_ptr<const qqrhdirect::Model> qqrhDirectModel;
};
struct Job {
    uint64_t generation = 0;
    int w = 0, h = 0;
    State state;
    fs::path output;
    Clock::time_point inputAt{};
    bool previewOnly = false, singlePass = false;
};
struct Box {
    int x, y, w, h;
    bool has(int px, int py) const { return px >= x && py >= y && px < x + w && py < y + h; }
    RECT rect() const { return {x, y, x + w, y + h}; }
};
struct Slider {
    int id;
    Box box;
    double lo, hi;
};
struct Options {
    bool cpu = false, test = false, gpuTest = false, compare = false, list = false, benchmark = false;
    bool engineProbe = false;
    engine_probe::Options probeOptions;
    rhgpu::EngineTuning gpuTuning;
    std::string interaction = "coarse";
    std::string previewDevice = "auto";
    std::string rhBackend = "auto";
    std::string renderer;
    std::string refine = "default";
    unsigned threads = 0;
    int width = 960, height = 540;
    fs::path render, view, uiSmoke, plotFive;
    Options() { gpuTuning.schedule=rhgpu::Schedule::Legacy; }
};
static unsigned rhCpuThreads(const Options& options) {
    return options.threads?options.threads:
        std::min(32u,std::max(1u,std::thread::hardware_concurrency()));
}
static bool rhPreferCpu(const Options& options,const State& state) {
    // Both favourite curvature views favoured CPU32 at every measured size.
    // Keep experimental GPU algebra opt-in until its target-machine results.
    return state.family==20 && (options.rhBackend=="cpu" ||
        (options.rhBackend=="auto" && state.generation==1));
}

class App {
    HWND window = nullptr;
    HFONT font = nullptr, smallFont = nullptr, titleFont = nullptr;
    float dpiScale = 1;
    int width = 1260, height = 946, canvasH = 526;
    bool help = false, dragging = false, orbitVisible = false, backValid = false, full = false;
    bool labControls = true;
    int focusBank = 1;
    int active = -1, dragX = 0, dragY = 0;
    Camera dragCamera;
    spin::Settings dragSpin;
    Box dragSpinPlot{0,0,1,1};
    State state, back;
    std::vector<catalogue::Entry> entries;
    std::vector<State> entryStates;
    bool inCatalogue = true, catalogueNotice = false;
    int selectedEntry = 0, activeEntry = 0, catalogueScroll = 0;
    Scene diagnostic;
    C picked{0, 0};
    C hilbertPicked{1, 0};
    bool hilbertInspectorVisible=true;
    std::string message = "Drag to pan. Wheel to zoom. Right-click to trace an orbit.";
    std::unique_ptr<Frame> frame;
    HDC paintDC = nullptr;
    HBITMAP paintBitmap = nullptr;
    HGDIOBJ paintPrevious = nullptr;
    int paintWidth = 0, paintHeight = 0;
    // UI-owned scheduling fields. Jobs copy their timing and quality before
    // entering the worker, which never reads the live interaction state.
    bool previewInFlight = false, deferredPreview = false, deferredColourOnly = true;
    bool settlePending = false, settleColourOnly = true, wheelInput = false;
    uint64_t previewGeneration = 0;
    Clock::time_point latestInputAt{}, lastSubmittedAt{}, settleDue{};
    std::optional<State> fullOrbitState;
    int fullOrbitWidth = 0, fullOrbitHeight = 0;
    std::mutex mutex;
    std::condition_variable cv;
    Job current, exportJob;
    bool pending = false, hqPending = false, closing = false, hqBusy = false;
    std::atomic<uint64_t> generation{0}, exportGeneration{0};
    std::thread worker;
    Options options;
    fs::path root;
    WINDOWPLACEMENT placement{};
    DWORD oldStyle = 0;
    std::string backend = "Starting Direct3D 11...";
    std::vector<Q> trailA, trailB;
    int sc(int n) const { return int(std::round(n * dpiScale)); }
    int cy(int n) const { return canvasH + sc(n); }
    bool spectralGen2() const { return state.family == 11 && state.generation == 2; }
    bool quantumMode() const { return state.family == 12; }
    bool sphereMode() const { return state.family == 13; }
    bool cubicMode() const { return state.family == 17; }
    bool quotientMode() const { return state.family == 18; }
    bool tripleMode() const { return state.family == 19; }
    bool masterMode() const { return state.family == 20; }
    bool extraMode() const { return state.family == 25; }
    bool hilbertMode() const { return extraMode() && state.generation>=9; }
    bool hilbertField() const { return extraMode() && frk_extra::is_field_driver(state.generation+17); }
    bool directQqrhMode() const { return state.family == 24; }
    bool connectedRhLambdaMode() const { return state.family == 22; }
    bool connectedQqrhLambdaMode() const { return state.family == 23; }
    bool connectedLambdaMode() const { return connectedRhLambdaMode() || connectedQqrhLambdaMode(); }
    bool legacyLambdaMode() const { return state.family == 21 || connectedLambdaMode(); }
    bool directAtlasMode() const { return cubicMode() || quotientMode() || tripleMode() || masterMode() || directQqrhMode() || extraMode(); }
    bool riemannMode() const { return state.family == 16; }
    bool projectMode() const { return state.family == 15; }
    bool spinMode() const { return state.family == 14; }
    bool spinAtlas() const { return spinMode() && state.generation == 3 && state.plane == 0; }
    bool sampleMode() const { return quantumMode() || sphereMode() || spinMode(); }
    int panelExtra() const { return spectralGen2() || sampleMode() || connectedLambdaMode() || directQqrhMode() || masterMode() ? 45 : 0; }
    Box legacyLambdaBox() const {
        const int side=std::max(1,std::min(width,canvasH));
        return {(width-side)/2,(canvasH-side)/2,side,side};
    }
    C legacyLambdaWorld(int x,int y) const {
        const auto b=legacyLambdaBox();
        return world(state.cameras[2],x-b.x,y-b.y,b.w,b.h);
    }
    Box row(int i, int count, int y, int h = 27) const {
        int gap = sc(5), margin = sc(10), w = (width - 2 * margin - (count - 1) * gap) / count;
        return {margin + i * (w + gap), cy(y), w, sc(h)};
    }
    std::vector<Slider> sliders() const {
        int m = sc(12), g = sc(10), w = (width - 2 * m - 2 * g) / 3;
        std::vector<Slider> result{{0, {m, cy(223), w, sc(38)}, -4, 4},
                                      {1, {m + w + g, cy(223), w, sc(38)}, -3, 3},
                                      {2, {m + 2 * (w + g), cy(223), w, sc(38)}, 32, 1024},
                                      {3, {m, cy(268), w, sc(38)}, -1, 3},
                                      {4, {m + w + g, cy(268), w, sc(38)}, -2, 2},
                                      {5, {m + 2 * (w + g), cy(268), w, sc(38)}, 0, 1},
                                      {6, {m, cy(313), w, sc(38)}, -4, 8},
                                      {7, {m + w + g, cy(313), w, sc(38)}, -4, 4},
                                      {8, {m + 2 * (w + g), cy(313), w, sc(38)}, .2, 3}};
        if (state.family == 3) {
            result[0].lo = result[1].lo = -.95;
            result[0].hi = result[1].hi = .95;
        }
        if (state.family == 6 || state.family == 8) {
            result[6].lo = -8; result[6].hi = 12;
            result[7].lo = -80; result[7].hi = 80;
        }
        if (state.family == 11) {
            result[0].lo = result[1].lo = -2; result[0].hi = result[1].hi = 2;
            result[3].lo = result[4].lo = -1; result[3].hi = result[4].hi = 1;
            result[6].lo = 0; result[6].hi = .49;
            result[7].lo = -2; result[7].hi = 2;
        }
        if (spectralGen2()) {
            result.push_back({9, {m, cy(358), w, sc(38)}, .01, .98});
            result.push_back({10, {m+w+g, cy(358), w, sc(38)}, 1, 3});
            result.push_back({11, {m+2*(w+g), cy(358), w, sc(38)}, -2, 2});
        }
        if (quantumMode()) {
            result[0].lo=0; result[0].hi=1;
            result[1].lo=-pi; result[1].hi=pi;
            result[2].lo=14; result[2].hi=23;
            result[3].lo=-pi; result[3].hi=pi;
            result[4].lo=-pi/2; result[4].hi=pi/2;
            result[6].lo=-pi; result[6].hi=pi;
            result[7].lo=-pi/2; result[7].hi=pi/2;
            result[8].lo=.25; result[8].hi=4;
            result.push_back({9,{m,cy(358),w,sc(38)},.4,8});
            result.push_back({10,{m+w+g,cy(358),w,sc(38)},0,3});
            result.push_back({11,{m+2*(w+g),cy(358),w,sc(38)},0,4294967295.0});
        }
        if (spinMode()) {
            result[0].lo=0; result[0].hi=20;
            result[1].lo=0; result[1].hi=pi;
            result[2].lo=spinAtlas() ? 5 : 14; result[2].hi=spinAtlas() ? 8 : 23;
            result[3].lo=0; result[3].hi=2;
            result[4].lo=6; result[4].hi=13;
            result[6].lo=-pi; result[6].hi=pi;
            result[7].lo=-pi/2; result[7].hi=pi/2;
            result[8].lo=.25; result[8].hi=4;
            result.push_back({9,{m,cy(358),w,sc(38)},.4,8});
            result.push_back({10,{m+w+g,cy(358),w,sc(38)},5,11});
            result.push_back({11,{m+2*(w+g),cy(358),w,sc(38)},0,4294967295.0});
        }
        if (sphereMode()) {
            result[0].lo=0; result[0].hi=.38;
            result[1].lo=-8; result[1].hi=0;
            result[2].lo=14; result[2].hi=23;
            result[3].lo=-pi; result[3].hi=pi;
            result[4].lo=-.85; result[4].hi=.85;
            result[5].lo=-pi; result[5].hi=pi;
            result[6].lo=-pi/2; result[6].hi=pi/2;
            result[7].lo=.25; result[7].hi=4;
            result[8].lo=.4; result[8].hi=8;
            result.push_back({9,{m,cy(358),w,sc(38)},0,1});
            result.push_back({10,{m+w+g,cy(358),w,sc(38)},0,1});
            result.push_back({11,{m+2*(w+g),cy(358),w,sc(38)},0,4294967295.0});
        }
        if(extraMode()) {
            result[0].lo=result[1].lo=-2.;result[0].hi=result[1].hi=2.;
            result[3].lo=-2.;result[3].hi=2.;
            if(hilbertMode()) {result[3].lo=.1;result[3].hi=1.25;}
        }
        if(directAtlasMode()) {
            result[0].lo=result[1].lo=-2;result[0].hi=result[1].hi=2;
        }
        if(masterMode() || directQqrhMode()) {
            result[3].lo=result[4].lo=-8;result[3].hi=result[4].hi=8;
            result[5].lo=-1;result[5].hi=1;
            result[6].lo=.025;result[6].hi=.49;
            result[7].lo=.16;result[7].hi=2.4;
            result[8].lo=3;result[8].hi=5;
        }
        if(masterMode()) {
            result[5].lo=.05;result[5].hi=3.;
            result.push_back({9,{m,cy(358),w,sc(38)},-12.,12.});
            if(labControls) {
                static const double lows[]={2,0,32,.05,0,.01,-pi,.1,.2};
                static const double highs[]={6,2,1024,1.95,.95,3,pi,3,8};
                for(int i=0;i<9;++i) {result[i].lo=lows[i];result[i].hi=highs[i];}
                result.push_back({10,{m+w+g,cy(358),w,sc(38)},-pi,pi});
                result.push_back({11,{m+2*(w+g),cy(358),w,sc(38)},0,1});
            }
        }
        if(directQqrhMode())result.push_back({9,{m,cy(358),w,sc(38)},state.generation==6 ? .25 : 0.,state.generation==6 ? 4. : 64.});
        if(connectedLambdaMode()) {
            result[0].lo=result[1].lo=-2;result[0].hi=result[1].hi=2;
            result[2].lo=32;result[2].hi=4096;
            result[3].lo=result[4].lo=-8;result[3].hi=result[4].hi=8;
            result[5].lo=-1;result[5].hi=1;
            result[6].lo=.025;result[6].hi=.49;
            result[7].lo=.16;result[7].hi=2.4;
            result[8].lo=3;result[8].hi=5;
            result.push_back({9,{m,cy(358),w,sc(38)},0,2});
            result.push_back({10,{m+w+g,cy(358),w,sc(38)},0,4});
            if(connectedQqrhLambdaMode())result.push_back({11,{m+2*(w+g),cy(358),w,sc(38)},0,64});
        } else if(legacyLambdaMode()) {
            result[0].lo=result[1].lo=-2;result[0].hi=result[1].hi=2;
            result[2].lo=32;result[2].hi=4096;
            result[3].lo=2;result[3].hi=8;
            result[4].lo=0;result[4].hi=4;
            result[5].lo=result[6].lo=result[7].lo=result[8].lo=0;
            result[5].hi=result[6].hi=result[7].hi=2;result[8].hi=4;
        }
        if(riemannMode()) {
            result[0].lo=1;result[0].hi=2;
            result[1].lo=-6;result[1].hi=-3;
            result[3].lo=result[4].lo=-2;result[3].hi=result[4].hi=2;
        }
        if (projectMode()) {
            if(state.generation==1) {
                result[3].lo=result[4].lo=-16; result[3].hi=result[4].hi=16;
                result[6].lo=12; result[6].hi=48;
                result[7].lo=2; result[7].hi=24;
            } else if(state.generation==2) {
                result[0].lo=16; result[0].hi=256;
                result[1].lo=0; result[1].hi=1024-state.project.piTerms;
            } else if(state.generation==3) {
                result[0].lo=4; result[0].hi=7;
                result[1].lo=0; result[1].hi=.01;
                result[3].lo=0; result[3].hi=3;
            } else if(state.generation==5) {
                result[0].lo=.02; result[0].hi=.9;
            }
        }
        return result;
    }
    bool enabled(int id) const {
        if(extraMode()) {
            const int driver=state.generation+17;
            if(hilbertField())return id==3;
            if(id==2)return true;
            if(id==3)return driver!=18 && driver!=21 && driver!=24;
            if(id==0 || id==1)return !frk_extra::seed_at_pixel(driver);
            return false;
        }
        if(directQqrhMode())return id>=0 && id<=9 && (id!=5 || (state.generation!=1 && state.generation!=6)) && (id!=9 || state.generation==5 || state.generation==6);
        if(connectedQqrhLambdaMode())return id>=0 && id<=11 && (id!=5 || state.generation!=1) && (id!=11 || state.generation==5);
        if(connectedRhLambdaMode())return id>=0 && id<=10 && (id!=5 || state.generation!=3);
        if(legacyLambdaMode()) return id>=0 && id<=8;
        if(masterMode()) return id>=0 && id<=(labControls ? 11 : 9) &&
            (labControls || state.plane!=2 || (id!=3 && id!=4));
        if(directAtlasMode()) return id==0 || id==1 || id==2 || id==5 || id==8;
        if(riemannMode()) return (state.generation<=3 && (id==0 || (id==1 && state.riemann.order==2))) ||
            id==2 || id==3 || id==4 || id==5 || id==8;
        if(projectMode()) {
            if(id==2 || id==5 || id==8) return true;
            if(state.generation==1) return id==0 || id==1 || id==3 || id==4 || id==6 || id==7;
            if(state.generation==2) return id==0 || id==1;
            if(state.generation==3) return id==0 || id==1 || id==3;
            return state.generation==5 && id==0;
        }
        if (spinMode()) {
            if (id==5 && !spinAtlas() && state.palette!=0) return false;
            return !spinAtlas() || (id!=6 && id!=7 && id!=9);
        }
        if (sphereMode()) {
            if (id==0) return state.generation==2;
            if (id==1 || id==3) return state.generation==3;
            return id!=10;
        }
        if (quantumMode()) return true;
        if (state.family == 9) return id == 5 || id == 8;
        if (state.family == 11) {
            if (id == 6) return true;
            if (id == 7 || id == 11) return state.plane != 0;
        }
        return !((state.plane == 1 && (id == 0 || id == 1)) || (state.plane == 2 && (id == 3 || id == 4)) ||
                 ((state.plane == 0 || (!state.family && state.orbitMode != 3)) && (id == 6 || id == 7)));
    }
    double value(int id) const {
        if(extraMode()) {
            if(id==0)return state.c.real();
            if(id==1)return state.c.imag();
            if(id==2)return state.iterations;
            if(id==3)return state.alpha.real();
            return 0;
        }
        if(directQqrhMode()) {
            switch(id) {
                case 0:return state.c.real();case 1:return state.c.imag();case 2:return state.iterations;
                case 3:return state.alpha.real();case 4:return state.alpha.imag();case 5:return state.seed.real();
                case 6:return state.spectralD;case 7:return state.seed.imag();case 8:return state.method+3;
                case 9:return state.generation==6 ? state.qqrhScale : state.qqrhSlope;default:return 0;
            }
        }
        if(connectedLambdaMode()) {
            switch(id) {
                case 0:return state.c.real();case 1:return state.c.imag();case 2:return state.iterations;
                case 3:return state.alpha.real();case 4:return state.alpha.imag();
                case 5:return connectedRhLambdaMode() && state.generation==3 ? -1. : state.seed.real();
                case 6:return state.spectralD;case 7:return state.seed.imag();case 8:return state.method+3;
                case 9:return state.bkqrQuality;case 10:return state.bkqrPrecision;
                case 11:return state.qqrhSlope;default:return 0;
            }
        }
        if(legacyLambdaMode()) {
            switch(id) {
                case 0:return state.c.real();case 1:return state.c.imag();case 2:return state.iterations;
                case 3:return state.alpha.real();case 4:return state.method;case 5:return state.orbitMode;
                case 6:return state.spectralOrder-1;case 7:return state.seed.real();case 8:return state.seed.imag();
                default:return 0;
            }
        }
        if(masterMode()) {
            if(labControls && id<9) {
                const auto &e=state.rhExperiment;
                switch(id) {
                    case 0:return e.power;case 1:return e.gain;case 2:return state.iterations;
                    case 3:return e.relaxation;case 4:return e.memory;case 5:return e.trapRadius;
                    case 6:return e.trapAngle;case 7:return e.exposure;case 8:return e.frequency;
                }
            }
            if(labControls && id==10)return state.rhExperiment.twist;
            if(labControls && id==11)return state.rhExperiment.blend;
            switch(id) {
                case 0:return state.c.real();case 1:return state.c.imag();case 2:return state.iterations;
                case 3:return state.alpha.real();case 4:return state.alpha.imag();
                case 5:return state.seed.real();
                case 6:return state.spectralD;case 7:return state.seed.imag();case 8:return state.spectralOrder+2;
                case 9:return std::isfinite(state.rhLevel) ? state.rhLevel : rhfresh::defaultSettings(state.generation-1).level;
                default:return 0;
            }
        }
        if(directAtlasMode()) {
            switch(id) {
                case 0:return state.c.real();case 1:return state.c.imag();case 2:return state.iterations;
                case 5:return state.hue;case 8:return state.bands;default:return 0;
            }
        }
        if(riemannMode()) {
            switch(id) {
                case 0:return state.riemann.order;case 1:return std::log10(state.riemann.anchor);
                case 2:return state.iterations;case 3:return state.alpha.real();case 4:return state.alpha.imag();
                case 5:return state.hue;case 8:return state.bands;default:return 0;
            }
        }
        if(projectMode()) {
            const auto &t=state.project;
            if(id==2) return state.iterations;
            if(id==5) return state.hue;
            if(id==8) return state.bands;
            if(state.generation==1) {
                if(id==0) return state.c.real(); if(id==1) return state.c.imag();
                if(id==3) return state.alpha.real(); if(id==4) return state.alpha.imag();
                if(id==6) return t.grhDegree; if(id==7) return t.grhRadius;
            }
            if(state.generation==2) return id==0 ? t.piTerms : id==1 ? t.piOffset : 0;
            if(state.generation==3) return id==0 ? t.nseModes : id==1 ? t.nseTime : id==3 ? t.nseSlope : 0;
            if(state.generation==5 && id==0) return t.pcdProbability;
            return 0;
        }
        if (spinMode()) {
            const auto &q=state.quantum;
            const auto &t=state.spin;
            switch(id) {
                case 0:return t.twist; case 1:return t.tilt;
                case 2:return spinAtlas() ? t.gridPower : q.samplePower;
                case 3:return t.damping; case 4:return t.burnPower; case 5:return state.hue;
                case 6:return q.yaw; case 7:return q.pitch; case 8:return q.exposure;
                case 9:return q.zoom; case 10:return t.windowPower; case 11:return q.randomSeed;
            }
        }
        if (sphereMode()) {
            const auto &q=state.quantum;
            switch(id) {
                case 0:return state.sphere.gap; case 1:return state.sphere.lambdaLog;
                case 2:return q.samplePower; case 3:return state.sphere.lambdaPhase;
                case 4:return state.sphere.boost; case 5:return q.yaw; case 6:return q.pitch;
                case 7:return q.exposure; case 8:return q.zoom; case 9:return state.hue;
                case 10:return 0; case 11:return q.randomSeed;
            }
        }
        if (quantumMode()) {
            const auto &q=state.quantum;
            switch(id) {
                case 0:return q.strength; case 1:return q.kick; case 2:return q.samplePower;
                case 3:return q.axisAzimuth; case 4:return q.axisElevation; case 5:return state.hue;
                case 6:return q.yaw; case 7:return q.pitch; case 8:return q.exposure;
                case 9:return q.zoom; case 10:return q.geometry; case 11:return q.randomSeed;
            }
        }
        switch (id) {
        case 0:
            return state.c.real();
        case 1:
            return state.c.imag();
        case 2:
            return state.iterations;
        case 3:
            return state.alpha.real();
        case 4:
            return state.alpha.imag();
        case 5:
            return state.hue;
        case 6:
            return state.family == 11 ? state.spectralD : state.seed.real();
        case 7:
            return state.family == 11 ? state.seed.real() : state.seed.imag();
        case 9:
            return state.spectralWeight;
        case 10:
            return state.spectralOrder;
        case 11:
            return state.seed.imag();
        default:
            return state.bands;
        }
    }
    std::string label(int id) const {
        if(directQqrhMode()) {
            if(id==5 && (state.generation==1 || state.generation==6))return "Angular x: unused by curvature";
            if(id==8)return "Theta Taylor degree = "+std::to_string(8<<state.method)+" [K]";
            if(id==9)return state.generation==6 ? "Scale kappa = "+num(state.qqrhScale,3) : state.generation==5 ? "Moving-boundary slope m = "+num(state.qqrhSlope,3) : "Slope m: moving-boundary channel only";
            static const char *names[]={"Selected c.real","Selected c.imag","Iteration budget","Fixed eta.real","Fixed eta.imag",
                "Angular x","Shift h","Radius scale b"};
            return std::string(names[id])+" = "+(id==2 ? std::to_string(state.iterations) : num(value(id),4));
        }
        if(connectedLambdaMode()) {
            if(id==5 && connectedQqrhLambdaMode() && state.generation==1)return "Angular x: unused by curvature";
            if(id==5 && connectedRhLambdaMode() && state.generation==3)return "Boundary x = -1 (fixed for Sub_2)";
            if(id==8)return "Theta Taylor degree = "+std::to_string(8<<state.method)+" [K]";
            if(id==9)return "Sampling quality = "+std::to_string(state.bkqrQuality);
            if(id==10)return std::string("Precision preset = ")+bkqrlambda::precisionName(state.bkqrPrecision);
            if(id==11)return state.generation==5 ? "Moving-boundary slope m = "+num(state.qqrhSlope,3) : "Slope m: moving-boundary channel only";
            static const char *names[]={"RH host c.real","RH host c.imag","Step budget","Inner eta.real","Inner eta.imag",
                "Angular x","Shift h","Radius scale b"};
            return std::string(names[id])+" = "+(id==2 ? std::to_string(state.iterations) : num(value(id),4));
        }
        if(legacyLambdaMode()) {
            if(id==4) return std::string("Host [F] = ")+bkqrlambda::hostName(state.method);
            if(id==5) return std::string("Branch [B] = ")+bkqrlambda::branchName(state.orbitMode);
            if(id==6) return std::string("Coordinates [C] = ")+bkqrlambda::couplingName(state.spectralOrder-1);
            if(id==7) return "Sampling quality = "+std::to_string(int(state.seed.real()));
            if(id==8) return std::string("Precision preset = ")+bkqrlambda::precisionName(int(state.seed.imag()));
            static const char *names[]={"Host c.real","Host c.imag","Step budget","Host power"};
            return std::string(names[id])+" = "+(id>=2 ? std::to_string(int(value(id))) : num(value(id),4));
        }
        if(masterMode()) {
            if(labControls) {
                static const char *names[]={"Host power","Feedback gain","Iteration budget","Relaxation","Memory",
                    "Trap radius","Trap angle (radians)","Exposure","Color frequency","Level",
                    "Twist (maps 9 / 11)","Blend (maps 10 / 11)"};
                return std::string(names[id])+" = "+((id==0 || id==2) ? std::to_string(int(value(id))) : num(value(id),4));
            }
            if(id==8) return "Theta Taylor degree = "+std::to_string(8<<(state.spectralOrder-1))+" [K]";
            if(state.plane==2 && id==3) return "lambda.real = pixel horizontal axis";
            if(state.plane==2 && id==4) return "lambda.imag = pixel vertical axis";
            static const char *names[]={"Selected c.real","Selected c.imag","Iteration budget","Feedback lambda.real",
                "Feedback lambda.imag","Normalized radius R","Shift h","Radius scale b","Taylor degree","Level"};
            return std::string(names[id])+" = "+(id==2 ? std::to_string(state.iterations) : num(value(id),4));
        }
        if(extraMode()) {
            if(id==3 && hilbertMode())return "Scaled radius R = "+num(state.alpha.real(),4)+" (r=200R)";
            if(hilbertField() && id==2)return "Direct Hilbert field / no iterations";
            if(id==0)return enabled(id)?"Selected c.real = "+num(state.c.real(),4):"c fixed / not used by root map";
            if(id==1)return enabled(id)?"Selected c.imag = "+num(state.c.imag(),4):"c fixed / not used by root map";
            if(id==2)return "Iteration budget = "+std::to_string(state.iterations);
            if(id==3)return enabled(id)?"Real deformation lambda = "+num(state.alpha.real(),4):"No deformation for this map";
            return "Original standalone controls / not used by this candidate";
        }
        if(tripleMode()) {
            if(id==3)return state.generation==1 ? "Conjugate cubic iteration" : state.generation==2 ? "Two alternating quadratic phases" : "Quadratic + moment pole";
            if(id==4)return "Coefficients fixed by RH moments";
            if(id==6)return "Atlas: all finite critical orbits";
            if(id==7)return "Finite nonescape: unresolved";
            static const char *names[]={"Selected c.real","Selected c.imag","Iteration budget","","","Palette shift","","","Band frequency"};
            if(id==2&&state.generation==2)return "Return cycles = "+std::to_string(state.iterations);
            return std::string(names[id])+" = "+(id==2 ? std::to_string(state.iterations) : num(value(id),4));
        }
        if(quotientMode()) {
            if(id==3) return "Direct rational: lambda*(z-b-1/z)+c";
            if(id==4) return "lambda and b fixed by RH moments";
            if(id==6) return "Critical +/-i; pole 0 -> infinity";
            if(id==7) return "Finite nonescape: unresolved";
            static const char *names[]={"Selected c.real", "Selected c.imag", "Iteration budget", "", "",
                "Palette shift", "", "", "Band frequency"};
            return std::string(names[id])+" = "+(id==2 ? std::to_string(state.iterations) : num(value(id),4));
        }
        if(cubicMode()) {
            if(id==3) return "Direct cubic: z^3 - A*z + c";
            if(id==4) return "A fixed by five RH moments";
            if(id==6) return "Atlas follows both critical points";
            if(id==7) return "Nonescape within budget: unresolved";
            static const char *names[]={"Selected c.real", "Selected c.imag", "Iteration budget", "", "",
                "Palette shift", "", "", "Band frequency"};
            return std::string(names[id])+" = "+(id==2 ? std::to_string(state.iterations) : num(value(id),4));
        }
        if(riemannMode()) {
            if(state.generation>=4 && id==0) return "H2 / five fixed moments";
            if(state.generation>=4 && id==1) return "Anchor: witness diagnostics only";
            if(state.generation>=4 && id==6) return state.generation==7 ? "Coordinate: Joukowski sheets" : "Coordinate: z=(q-mu)/sigma";
            if(state.generation>=4 && id==7) return "Fixed moment polynomial / image metadata";
            if(id==0) return "Finite moment order = H"+std::to_string(state.riemann.order)+
                (state.riemann.order==1 ? " / 2 nodes" : " / 3 nodes");
            if(id==1) return state.riemann.order==2 ? "log10 anchor = "+num(value(id),4) : "Anchor: H2 only";
            if(id==6) return state.generation==1 ? "Coordinate: q / q_max" : "Cubic: (27/4) z (1-z)^2";
            if(id==7) return "Witness data saved with each image";
            static const char *names[]={"Order", "Anchor", "Iteration budget", "Relaxation.real", "Relaxation.imag",
                "Palette shift", "Coordinate", "Witness", "Band frequency"};
            return std::string(names[id])+" = "+(id==2 ? std::to_string(state.iterations) : num(value(id),4));
        }
        if(projectMode()) {
            if(!enabled(id)) return "Fixed by this project map";
            std::array<std::string,9> names={"Offset a.real", "Offset a.imag", "Iteration budget", "Multiplier b.real",
                "Multiplier b.imag", "Palette shift", "GRH approximation degree", "GRH coefficient scale", "Band frequency"};
            if(state.generation==2) {names[0]="Actual pi digit pairs";names[1]="Pair offset from 1024";}
            if(state.generation==3) {names[0]="NSE mode count";names[1]="Diffusion time nu*t";names[3]="Energy slope";}
            if(state.generation==5) names[0]="q=2 probability p";
            const bool integral=id==2 || (state.generation==1 && id==6) || (state.generation==2 && id<2) || (state.generation==3 && id==0);
            return names[id]+" = "+(integral ? std::to_string(int(value(id))) : num(value(id),state.generation==3 && id==1 ? 6 : 4));
        }
        if (spinMode()) {
            static const char *names[]={"Twist k", "Tilt beta", "Sample budget", "Damping tau", "Burn-in",
                "Palette shift", "View yaw", "View pitch", "Exposure", "View zoom", "Lyapunov window", "Initial seed"};
            if (id==0 || id==3) return std::string(spinAtlas() ? "Selected " : "")+names[id]+" = "+num(value(id),4);
            if (id==2) return std::string(spinAtlas() ? "Grid width" : "Sample budget")+" = 2^"+
                std::to_string(int(value(id)))+" / "+std::to_string(uint64_t(1)<<int(value(id)));
            if (id==4 || id==10) return std::string(names[id])+" = 2^"+std::to_string(int(value(id)))+
                " / "+std::to_string(uint64_t(1)<<int(value(id)));
            if (id==11) return std::string(names[id])+" = "+std::to_string(state.quantum.randomSeed);
            if (!enabled(id)) return std::string(names[id])+(id==5 ? ": fixed signed colours" : ": sphere views");
            return std::string(names[id])+" = "+num(value(id),3);
        }
        if (sphereMode()) {
            static const char *names[]={"Circle separation", "log10 |lambda|", "Sample budget",
                "Lambda phase (radians)", "Conformal display boost", "View yaw", "View pitch",
                "Density exposure", "View zoom", "Palette shift", "Construction", "Replay seed"};
            if (id==10) return state.generation==3 ? "Construction: 6 inverse branches" : "Construction: tetrahedral circles";
            if (!enabled(id)) return std::string(names[id])+(id==0 ? ": Schottky only" : ": McMullen only");
            if (id==2) return std::string(names[id])+" = 2^"+std::to_string(state.quantum.samplePower)+
                " / "+std::to_string(uint64_t(1)<<state.quantum.samplePower);
            if (id==11) return std::string(names[id])+" = "+std::to_string(state.quantum.randomSeed);
            return std::string(names[id])+" = "+num(value(id),3);
        }
        if (quantumMode()) {
            static const char *names[]={"Measurement strength", "Unitary kick (radians)", "Sample budget",
                "Kick axis azimuth", "Kick axis elevation", "Palette shift", "View yaw", "View pitch",
                "Density exposure", "View zoom", "Measurement geometry", "Replay seed"};
            if (id==2) return std::string(names[id])+" = 2^"+std::to_string(state.quantum.samplePower)+
                " / "+std::to_string(uint64_t(1)<<state.quantum.samplePower);
            if (id==10) return std::string(names[id])+" = "+quantum::geometryName(state.quantum.geometry);
            if (id==11) return std::string(names[id])+" = "+std::to_string(state.quantum.randomSeed);
            return std::string(names[id])+" = "+num(value(id),3);
        }
        std::array<std::string, 12> names = {"Offset c.real",   "Offset c.imag",   "Iteration budget",
                                            "Multiplier.real", "Multiplier.imag", "Palette shift",
                                            "Seed.real",       "Seed.imag",       "Band frequency",
                                            "Middle atom weight", "Jet order s", "Parameter seed.imag"};
        if (state.family == 1) {
            names[0] = "Branch a.real";
            names[1] = "Branch a.imag";
            names[3] = "Branch b.real";
            names[4] = "Branch b.imag";
        }
        if (state.family == 4) {
            names[0] = "Feedback.real";
            names[1] = "Feedback.imag";
            names[3] = "Relaxation.real";
            names[4] = "Relaxation.imag";
        }
        if (state.family == 3) {
            names[0] = "Embedding q.real";
            names[1] = "Embedding q.imag";
            names[3] = "Relaxation.real";
            names[4] = "Relaxation.imag";
        }
        if (state.family == 2 || state.family >= 5) {
            names[0] = "Target.real";
            names[1] = "Target.imag";
            names[3] = "Relaxation.real";
            names[4] = "Relaxation.imag";
        }
        if (state.family == 7 || state.family == 8) {
            names[0] = "Normalized target.real";
            names[1] = "Normalized target.imag";
        }
        if (state.family == 9) {
            names[0] = "Fixed by observable"; names[1] = "Fixed by observable";
            names[2] = "Direct diagnostic"; names[3] = "Fixed by observable"; names[4] = "Fixed by observable";
            names[5] = "Phase palette shift"; names[6] = "Shared parameter plane"; names[7] = "Shared parameter plane";
            names[8] = "Contour / phase frequency";
        }
        if (state.family == 10) {
            names[0] = "Phase kick kappa"; names[1] = "Phase twist tau";
            names[3] = "Newton relaxation.real"; names[4] = "Newton relaxation.imag";
            names[6] = "Seed.real"; names[7] = "Seed.imag"; names[8] = "Flow ribbon frequency";
        }
        if (state.family == 11) {
            names[0] = "Offset c.real"; names[1] = "Offset c.imag";
            names[3] = "Feedback eta.real"; names[4] = "Feedback eta.imag";
            names[6] = "Spectral spacing d"; names[7] = "Parameter seed.real";
            names[8] = "Colour frequency";
        }
        if (!enabled(id)) {
            if (id == 6 || id == 7 || id == 11)
                return names[id] +
                       (state.family ? ": used in parameter planes" : ": select Seed in parameter plane");
            return names[id] + ": pixel " + ((id == 0 || id == 3) ? "horizontal" : "vertical") + " axis";
        }
        return names[id] + " = " + (id == 2 ? std::to_string(state.iterations) : id == 10 ? std::to_string(state.spectralOrder) : num(value(id), 3));
    }
    void setValue(int id, double v) {
        if (!enabled(id))
            return;
        auto sl = sliders()[id];
        v = std::clamp(v, sl.lo, sl.hi);
        if(extraMode()) {
            if(id==0)state.c.real(v);
            else if(id==1)state.c.imag(v);
            else if(id==2)state.iterations=std::clamp(int(std::round(v/16))*16,32,1024);
            else if(id==3)state.alpha.real(v);
            state.preset=-1;refresh();request();return;
        }
        if(directQqrhMode()) {
            if(id==0)state.c.real(v);else if(id==1)state.c.imag(v);
            else if(id==2)state.iterations=std::clamp(int(std::round(v/16))*16,32,1024);
            else if(id==3)state.alpha.real(v);else if(id==4)state.alpha.imag(v);
            else if(id==5)state.seed.real(v);else if(id==6)state.spectralD=v;else if(id==7)state.seed.imag(v);
            else if(id==8)state.method=std::clamp(int(std::round(v))-3,0,2);
            else if(id==9) {
                if(state.generation==6)state.qqrhScale=std::clamp(v,.25,4.);
                else state.qqrhSlope=std::clamp(v,0.,64.);
            }
            state.preset=-1;refresh();request();return;
        }
        if(connectedLambdaMode()) {
            if(id==0)state.c.real(v);else if(id==1)state.c.imag(v);
            else if(id==2)state.iterations=std::clamp(int(std::round(v/16))*16,32,4096);
            else if(id==3)state.alpha.real(v);else if(id==4)state.alpha.imag(v);
            else if(id==5)state.seed.real(v);else if(id==6)state.spectralD=v;else if(id==7)state.seed.imag(v);
            else if(id==8)state.method=std::clamp(int(std::round(v))-3,0,2);
            else if(id==9)state.bkqrQuality=std::clamp(int(std::round(v)),0,2);
            else if(id==10)state.bkqrPrecision=std::clamp(int(std::round(v)),0,4);
            else if(id==11)state.qqrhSlope=std::clamp(v,0.,64.);
            state.preset=-1;orbitVisible=false;refresh();request();return;
        }
        if(legacyLambdaMode()) {
            if(id==0)state.c.real(v);else if(id==1)state.c.imag(v);
            else if(id==2)state.iterations=std::clamp(int(std::round(v/16))*16,32,4096);
            else if(id==3)state.alpha=C(std::clamp(int(std::round(v)),2,8),0);
            else if(id==4)state.method=std::clamp(int(std::round(v)),0,4);
            else if(id==5)state.orbitMode=std::clamp(int(std::round(v)),0,2);
            else if(id==6)state.spectralOrder=std::clamp(int(std::round(v)),0,2)+1;
            else if(id==7)state.seed.real(std::clamp(int(std::round(v)),0,2));
            else if(id==8)state.seed.imag(std::clamp(int(std::round(v)),0,4));
            state.preset=-1;orbitVisible=false;refresh();request();return;
        }
        if(masterMode()) {
            if(labControls && (id==10 || id==11)) {
                if(id==10)state.rhExperiment.twist=v;else state.rhExperiment.blend=v;
                state.preset=-1;refresh();request();return;
            }
            if(labControls && id<9) {
                auto &e=state.rhExperiment;
                switch(id) {
                    case 0:e.power=std::clamp(int(std::round(v)),2,6);break;
                    case 1:e.gain=v;break;
                    case 2:state.iterations=std::clamp(int(std::round(v/16))*16,32,1024);break;
                    case 3:e.relaxation=v;break;case 4:e.memory=v;break;case 5:e.trapRadius=v;break;
                    case 6:e.trapAngle=v;break;case 7:e.exposure=v;break;case 8:e.frequency=v;break;
                }
                state.preset=-1;refresh();request(id==7 || id==8);return;
            }
            if(id==0) state.c.real(v);
            else if(id==1) state.c.imag(v);
            else if(id==2) state.iterations=std::clamp(int(std::round(v/16))*16,32,1024);
            else if(id==3) state.alpha.real(v);
            else if(id==4) state.alpha.imag(v);
            else if(id==5) state.seed.real(v);
            else if(id==6) state.spectralD=v;
            else if(id==7) state.seed.imag(v);
            else if(id==8) state.spectralOrder=std::clamp(int(std::round(v))-2,1,3);
            else if(id==9) state.rhLevel=v;
            if(id!=2) state.preset=-1;
            refresh();request();return;
        }
        if(directAtlasMode()) {
            if(id==0) state.c.real(v);
            else if(id==1) state.c.imag(v);
            else if(id==2) state.iterations=std::clamp(int(std::round(v/16))*16,32,1024);
            else if(id==5) state.hue=v;
            else if(id==8) state.bands=v;
            if(id<2) state.preset=-1;
            refresh();request();return;
        }
        if(riemannMode()) {
            if(id==0) state.riemann.order=std::clamp(int(std::round(v)),1,2);
            else if(id==1) state.riemann.anchor=std::clamp(std::pow(10.,v),1e-6,1e-3);
            else if(id==2) state.iterations=std::clamp(int(std::round(v/16))*16,32,1024);
            else if(id==3 || id==4) {
                C relaxation=state.alpha;
                if(id==3) relaxation.real(v);else relaxation.imag(v);
                const double magnitude=std::abs(relaxation);
                if(magnitude>2) relaxation*=(2-1e-12)/magnitude;
                if(std::abs(relaxation)<1e-12) relaxation=C(.001,0);
                state.alpha=relaxation;
            } else if(id==5) state.hue=v;
            else if(id==8) state.bands=v;
            if(id!=2 && id!=5 && id!=8) state.preset=-1;
            refresh();request();return;
        }
        if(projectMode()) {
            auto &t=state.project;
            if(id==2) state.iterations=std::clamp(int(std::round(v/16))*16,32,1024);
            else if(id==5) state.hue=v;
            else if(id==8) state.bands=v;
            else if(state.generation==1) {
                if(id==0) state.c.real(v); else if(id==1) state.c.imag(v);
                else if(id==3) state.alpha.real(v); else if(id==4) state.alpha.imag(v);
                else if(id==6) t.grhDegree=std::clamp(int(std::round(v/2))*2,12,48);
                else if(id==7) t.grhRadius=v;
                if(std::abs(state.alpha)<1e-12) state.alpha=C(.001,0);
            } else if(state.generation==2) {
                if(id==0) {t.piTerms=std::clamp(int(std::round(v)),16,256);t.piOffset=std::min(t.piOffset,1024-t.piTerms);}
                else if(id==1) t.piOffset=std::clamp(int(std::round(v)),0,1024-t.piTerms);
            } else if(state.generation==3) {
                if(id==0) t.nseModes=std::clamp(int(std::round(v)),4,7);
                else if(id==1) t.nseTime=v; else if(id==3) t.nseSlope=v;
            } else if(state.generation==5 && id==0) t.pcdProbability=v;
            if(id!=2 && id!=5 && id!=8) state.preset=-1;
            refresh();request();return;
        }
        if (spinMode()) {
            auto &q=state.quantum;
            auto &t=state.spin;
            switch(id) {
                case 0:t.twist=v;break; case 1:t.tilt=v;break;
                case 2:if(spinAtlas()) t.gridPower=std::clamp(int(std::round(v)),5,8);
                       else q.samplePower=std::clamp(int(std::round(v)),14,23);break;
                case 3:t.damping=v;break; case 4:t.burnPower=std::clamp(int(std::round(v)),6,13);break;
                case 5:state.hue=v;break; case 6:q.yaw=v;break; case 7:q.pitch=v;break;
                case 8:q.exposure=v;break; case 9:q.zoom=v;break;
                case 10:t.windowPower=std::clamp(int(std::round(v)),5,11);break;
                case 11:q.randomSeed=uint32_t(std::round(v));break;
            }
            if (id<5 || id==10 || id==11) state.preset=-1;
            refresh();request();return;
        }
        if (sphereMode()) {
            auto &q=state.quantum;
            switch(id) {
                case 0:state.sphere.gap=v;break; case 1:state.sphere.lambdaLog=v;break;
                case 2:q.samplePower=std::clamp(int(std::round(v)),14,23);break;
                case 3:state.sphere.lambdaPhase=v;break; case 4:state.sphere.boost=v;break;
                case 5:q.yaw=v;break; case 6:q.pitch=v;break; case 7:q.exposure=v;break;
                case 8:q.zoom=v;break; case 9:state.hue=v;break;
                case 11:q.randomSeed=uint32_t(std::round(v));break;
            }
            if (id==0 || id==1 || id==3 || id==4 || id==11) state.preset=-1;
            refresh();request();return;
        }
        if (quantumMode()) {
            auto &q=state.quantum;
            switch(id) {
                case 0:q.strength=v;break; case 1:q.kick=v;break;
                case 2:q.samplePower=std::clamp(int(std::round(v)),14,23);break;
                case 3:q.axisAzimuth=v;break; case 4:q.axisElevation=v;break;
                case 5:state.hue=v;break; case 6:q.yaw=v;break; case 7:q.pitch=v;break;
                case 8:q.exposure=v;break; case 9:q.zoom=v;break;
                case 10:q.geometry=std::clamp(int(std::round(v)),0,3);break;
                case 11:q.randomSeed=uint32_t(std::round(v));break;
            }
            if (id<5 || id==10 || id==11) state.preset=-1;
            refresh(); request(); return;
        }
        switch (id) {
        case 0:
            state.c.real(v);
            break;
        case 1:
            state.c.imag(v);
            break;
        case 2:
            state.iterations = std::clamp(int(std::round(v / 16)) * 16, 32, 1024);
            break;
        case 3:
            state.alpha.real(v);
            break;
        case 4:
            state.alpha.imag(v);
            break;
        case 5:
            state.hue = v;
            break;
        case 6:
            if (state.family == 11) state.spectralD = v;
            else state.seed.real(v);
            break;
        case 7:
            if (state.family == 11) state.seed.real(v);
            else state.seed.imag(v);
            break;
        case 8:
            state.bands = v;
            break;
        case 9:
            state.spectralWeight = v;
            break;
        case 10:
            state.spectralOrder = std::clamp(int(std::round(v)), 1, 3);
            break;
        case 11:
            state.seed.imag(v);
            break;
        }
        if ((id < 5 && id != 2) || (state.family == 11 && (id == 6 || id == 7 || id >= 9)))
            state.preset = -1;
        refresh();
        request();
    }
    void moveSlider(int x) {
        auto sl = sliders()[active];
        double t = double(x - sl.box.x - sc(8)) / std::max(20, sl.box.w - sc(16));
        setValue(active, sl.lo + std::clamp(t, 0., 1.) * (sl.hi - sl.lo));
    }
    void notice(std::string text) {
        catalogueNotice =
            text.rfind("Error:", 0) == 0 || text.rfind("Render error:", 0) == 0 || text.rfind("HQ ", 0) == 0;
        message = std::move(text);
        InvalidateRect(window, nullptr, FALSE);
    }
    void refresh() {
        canvasH = std::max(sc(120), height - sc(422 + panelExtra()));
        diagnostic = scene(state);
        trailA.clear();
        trailB.clear();
        if (!orbitVisible || state.family == 9 || sampleMode() || legacyLambdaMode())
            return;
        if(extraMode()) {
            if(state.plane!=1 || hilbertField())return;
            frk_extra::Params q;q.lambda=state.alpha.real();q.julia_c=state.c;
            q.julia=true;q.iterations=state.iterations;
            auto orbit=frk_extra::start(state.generation+17,picked,q);
            for(int j=0;j<std::min(state.iterations,256);++j) {
                if(!frk_extra::finite(orbit.z) || std::abs(orbit.z)>1e9)break;
                trailA.push_back(point(orbit.z));
                orbit=frk_extra::step(state.generation+17,orbit,q.lambda);
            }
            return;
        }
        if(directQqrhMode()) {
            if(state.plane!=1 || !diagnostic.qqrhDirectModel || !diagnostic.qqrhDirectModel->validModel)return;
            std::vector<C> trail;
            diagnostic.qqrhDirectModel->orbit(picked,state.c,state.iterations,&trail);
            for(C z:trail) {
                if(!std::isfinite(z.real()) || !std::isfinite(z.imag()) || std::abs(z)>=1e100)break;
                trailA.push_back(point(z));
            }
            return;
        }
        if(masterMode()) {
            if(state.plane!=1 || !diagnostic.masterModel || !diagnostic.masterModel->validModel) return;
            std::vector<C> trail;
            diagnostic.masterModel->orbitAtLambda(masterPhysicalPicked(),state.c,state.alpha,state.iterations,&trail);
            C previous=picked;
            for(C z:trail) {
                if(!std::isfinite(z.real()) || !std::isfinite(z.imag()) || std::abs(z)>=1e100) break;
                C chartPoint=rhmaster::chartInverse(z,state.rhExperiment.chart,previous);
                if(!std::isfinite(chartPoint.real()) || !std::isfinite(chartPoint.imag()) || std::abs(chartPoint)>=1e100)break;
                trailA.push_back(point(chartPoint));previous=chartPoint;
            }
            return;
        }
        if(tripleMode()) {
            if(state.plane!=1 || !diagnostic.tripleModel || !diagnostic.tripleModel->validModel) return;
            std::vector<C> trail;
            diagnostic.tripleModel->orbit(picked,state.c,state.iterations,&trail);
            for(C z:trail) {
                if(!std::isfinite(z.real()) || !std::isfinite(z.imag()) || std::abs(z)>=1e100) break;
                trailA.push_back(point(z));
            }
            return;
        }
        if(quotientMode()) {
            if(state.plane!=1 || !diagnostic.quotientModel || !diagnostic.quotientModel->validModel) return;
            std::vector<C> trail;
            diagnostic.quotientModel->orbit(picked,state.c,state.iterations,&trail);
            for(C z:trail) {
                if(!std::isfinite(z.real()) || !std::isfinite(z.imag()) || std::abs(z)>=1e100) break;
                trailA.push_back(point(z));
            }
            return;
        }
        if(cubicMode()) {
            if(state.plane!=1 || !diagnostic.cubicModel || !diagnostic.cubicModel->validModel) return;
            std::vector<C> trail;
            diagnostic.cubicModel->orbit(picked,state.c,state.iterations,&trail);
            for(C z:trail) {
                if(!std::isfinite(z.real()) || !std::isfinite(z.imag()) || std::abs(z)>=1e100) break;
                trailA.push_back(point(z));
            }
            return;
        }
        if(riemannMode()) {
            if(!diagnostic.riemannModel || !diagnostic.riemannModel->validModel) return;
            std::vector<C> trail;
            diagnostic.riemannModel->orbit(picked,std::min(state.iterations,192),&trail);
            for(C z:trail) if(std::isfinite(z.real()) && std::isfinite(z.imag()) && std::abs(z)<1e100)
                trailA.push_back(point(z));
            return;
        }
        if(projectMode()) {
            std::vector<C> trail;
            diagnostic.projectModel->orbit(picked,std::min(state.iterations,192),&trail);
            for(C z:trail) if(std::isfinite(z.real()) && std::isfinite(z.imag()) && std::abs(z)<1e100)
                trailA.push_back(point(z));
            return;
        }
        if (state.family == 11) {
            std::vector<C> trail;
            if (spectralGen2())
                spectral::orbitGen2(state.plane == 0 ? picked : state.seed, state.c, state.alpha,
                                   state.spectralD, state.spectralWeight, state.spectralOrder,
                                   std::min(state.iterations, 192), &trail);
            else
                spectral::orbit(state.plane == 0 ? picked : state.seed, state.c, state.alpha,
                                state.spectralD, std::min(state.iterations, 192), &trail);
            for (C z : trail)
                if (std::isfinite(z.real()) && std::isfinite(z.imag()) && std::abs(z) < 1e100)
                    trailA.push_back(point(z));
            return;
        }
        if (state.family) {
            std::vector<C> trail;
            garden::orbit(state.plane == 0 ? picked : state.seed, state.c, state.alpha,
                          *diagnostic.gardenData, state.method, std::min(state.iterations, 192), &trail);
            for (C z : trail)
                if (std::isfinite(z.real()) && std::isfinite(z.imag()))
                    trailA.push_back(point(z));
            return;
        }
        auto add = [&](std::vector<Q> &list, Q q) {
            for (int i = 0; i < 96; i++) {
                list.push_back(q);
                q = advance(q, state.c, state.alpha);
            }
        };
        if (state.plane == 0)
            add(trailA, point(picked));
        else if (state.orbitMode == 3)
            add(trailA, point(state.seed));
        else {
            if (state.orbitMode != 2)
                add(trailA, point({2.4, .8}));
            if (state.orbitMode != 1)
                add(trailB, point({2.4, -.8}));
        }
    }
    static bool sameOrbitScene(const State& a, const State& b) {
        if(a.family!=20 || b.family!=20 || a.generation!=b.generation || a.plane!=b.plane ||
           a.iterations!=b.iterations || a.spectralOrder!=b.spectralOrder || a.method!=b.method ||
           a.c!=b.c || a.alpha!=b.alpha || a.seed!=b.seed || a.spectralD!=b.spectralD ||
           !(a.rhLevel==b.rhLevel || (std::isnan(a.rhLevel)&&std::isnan(b.rhLevel)))) return false;
        for(int i=0;i<3;++i)
            if(a.cameras[i].x!=b.cameras[i].x || a.cameras[i].y!=b.cameras[i].y ||
               a.cameras[i].span!=b.cameras[i].span) return false;
        auto left=a.rhExperiment, right=b.rhExperiment;
        // These controls affect the final paint only. Keep chart, metric and
        // trap geometry in the key because they can require different data.
        left.exposure=right.exposure; left.frequency=right.frequency;
        return rhmaster::sameExperiment(left,right);
    }
    static bool samePaintScene(const State& a, const State& b) {
        return sameOrbitScene(a,b) && a.palette==b.palette && a.hue==b.hue && a.bands==b.bands &&
            a.rhExperiment.exposure==b.rhExperiment.exposure &&
            a.rhExperiment.frequency==b.rhExperiment.frequency;
    }
    void submitRender(bool previewOnly, bool colourOnly, Clock::time_point inputAt,
                      bool finalPass=false) {
        if (!window || canvasH < 1 || width < 1)
            return;
        Job j;
        j.generation = ++generation;
        j.w = width;
        j.h = canvasH;
        if (inCatalogue) {
            auto box = previewBox();
            j.w = std::min(box.w, 384);
            j.h = std::max(32, j.w * box.h / box.w);
        }
        j.state = state;
        j.inputAt = inputAt;
        j.previewOnly = previewOnly;
        // An explicit paint hint is never enough by itself: a stale full
        // frame must not turn a newer dynamics change into a recolour job.
        j.singlePass = finalPass || (colourOnly && fullOrbitState &&
            fullOrbitWidth==j.w && fullOrbitHeight==j.h && sameOrbitScene(*fullOrbitState,j.state));
        lastSubmittedAt=Clock::now();
        if(previewOnly) {
            previewInFlight=true;
            previewGeneration=j.generation;
        }
        {
            std::lock_guard<std::mutex> lock(mutex);
            current = j;
            pending = true;
        }
        cv.notify_one();
        InvalidateRect(window, nullptr, FALSE);
    }
    void interactionTick() {
        if(!masterMode() || inCatalogue) {
            deferredPreview=settlePending=previewInFlight=false;
            KillTimer(window,INTERACTION_TIMER);
            return;
        }
        const auto now=Clock::now();
        // Let one preview finish even under a stream of mouse events. New
        // input replaces only the pending snapshot instead of cancelling
        // every in-flight frame before it can become visible.
        if(deferredPreview && !previewInFlight &&
           now-lastSubmittedAt>=std::chrono::milliseconds(33)) {
            const bool colourOnly=deferredColourOnly;
            deferredPreview=false; deferredColourOnly=true;
            submitRender(true,colourOnly,latestInputAt);
        }
        if(settlePending && !deferredPreview && !previewInFlight &&
           !dragging && active<0 && now>=settleDue) {
            settlePending=false;
            // A compatible full recolour may already have completed while
            // the slider was held. Do not repeat that exact paint on release.
            if(!frame || frame->stride!=1 || frame->targetW!=width ||
               frame->targetH!=canvasH || !samePaintScene(frame->state,state))
                submitRender(false,settleColourOnly,latestInputAt,true);
        }
        if(!deferredPreview && !settlePending && !previewInFlight)
            KillTimer(window,INTERACTION_TIMER);
    }
    void finishInteraction() {
        if(!masterMode() || inCatalogue || (!settlePending && !deferredPreview)) return;
        settleDue=Clock::now()+std::chrono::milliseconds(100);
        SetTimer(window,INTERACTION_TIMER,16,nullptr);
        interactionTick();
    }
    void request(bool colourOnly=false) {
        if(!window || canvasH<1 || width<1) return;
        latestInputAt=Clock::now();
        if(masterMode() && !inCatalogue && (dragging || active>=0 || wheelInput)) {
            if(!settlePending)settleColourOnly=colourOnly;
            else settleColourOnly=settleColourOnly&&colourOnly;
            settlePending=true;
            settleDue=latestInputAt+std::chrono::milliseconds(100);
            if(!deferredPreview)deferredColourOnly=colourOnly;
            else deferredColourOnly=deferredColourOnly&&colourOnly;
            deferredPreview=true;
            SetTimer(window,INTERACTION_TIMER,16,nullptr);
            if(options.interaction=="latest") {
                deferredPreview=false; deferredColourOnly=true;
                submitRender(true,colourOnly,latestInputAt);
                InvalidateRect(window,nullptr,FALSE);
                return;
            }
            interactionTick();
            InvalidateRect(window,nullptr,FALSE);
            return;
        }
        // Discrete commands, view loads and catalogue changes replace all
        // queued interaction work and cancel the old revision immediately.
        deferredPreview=settlePending=previewInFlight=false;
        deferredColourOnly=settleColourOnly=true;
        KillTimer(window,INTERACTION_TIMER);
        submitRender(false,colourOnly,latestInputAt);
    }
    void workerLoop() {
        Renderer cpuPreviewEngine;
        // The measured small RH frames favour CPU32 over GPU16 plus repairs.
        // This renderer owns only a lazy CPU pool; it creates no D3D device.
        const unsigned previewThreads=rhCpuThreads(options);
        const unsigned compilingThreads=std::max(1u,std::min(8u,previewThreads/2));
        cpuPreviewEngine.init(true,previewThreads);
        rhbackground::BackgroundGpu gpuBackend;
        gpuBackend.start(options.cpu,options.threads,[handle=window,requestedCpu=options.cpu](bool ok) {
            if(!requestedCpu)PostMessageW(handle,GPU_READY,ok?1:0,0);
        });
        bool compilerReported=options.cpu;
        std::cout << "BACKEND=" << (options.cpu?"CPU requested":"CPU while GPU shaders compile in background")
                  << "\nDETAIL=UI logic and CPU rendering continue independently of shader compilation"
                  << "\nCPU_THREADS=" << (options.threads ? std::to_string(options.threads) : "auto")
                  << "\nREFINE_PROFILE=" << options.refine
                  << "\nINTERACTION_POLICY=" << options.interaction
                  << "\nPREVIEW_DEVICE=" << options.previewDevice
                  << "\nPREVIEW_CPU_THREADS=" << previewThreads
                  << "\nCPU_WORKERS_DURING_COMPILE=" << compilingThreads
                  << "\nRH_BACKEND_POLICY=" << options.rhBackend
                  << "\nGPU_WAIT=" << (options.gpuTuning.waitPolicy==rhgpu::WaitPolicy::Responsive?"responsive":"legacy-sleep")
                  << "\nGPU_SCHEDULE=" << (options.gpuTuning.schedule==rhgpu::Schedule::Legacy?"legacy":"adaptive") << std::endl;
        auto startup = new std::string("Ready. The active family backend is shown below the canvas. H opens controls.");
        PostMessageW(window, NOTICE, 0, reinterpret_cast<LPARAM>(startup));
        for (;;) {
            Job job;
            bool exporting = false;
            {
                std::unique_lock<std::mutex> lock(mutex);
                cv.wait(lock, [&] { return closing || pending || hqPending; });
                if (closing)
                    break;
                if (hqPending) {
                    job = exportJob;
                    hqPending = false;
                    exporting = true;
                } else {
                    job = current;
                    pending = false;
                }
            }
            try {
                Renderer* accelerated=gpuBackend.renderer();
                if(accelerated)accelerated->setGpuTuning(options.gpuTuning);
                cpuPreviewEngine.setThreadLimit(gpuBackend.loading()?compilingThreads:previewThreads);
                if(!compilerReported&&gpuBackend.ready()) {
                    compilerReported=true;
                    std::cout<<"BACKGROUND_GPU="<<(accelerated?"READY":"FAILED")
                        <<" COMPILE_MS="<<gpuBackend.elapsedMs()
                        <<" DETAIL="<<(accelerated?accelerated->backend:gpuBackend.error())<<std::endl;
                }
                Scene r = scene(job.state);
                if (exporting) {
                    Renderer& exportEngine=(!accelerated||rhPreferCpu(options,r.state))?cpuPreviewEngine:*accelerated;
                    std::vector<uint32_t> data;
                    bool ok = exportEngine.render(r, job.w, job.h, double(job.h) / job.w, 2, data,
                                            &exportGeneration, job.generation);
                    if (ok) {
                        writeBmp(job.output, data, job.w, job.h);
                        const auto qs=exportEngine.quantumStats();
                        const auto ss=exportEngine.sphericalStats();
                        const auto ts=exportEngine.spinStats();
                        metadata(job.output, job.state, job.w, job.h, 2, exportEngine.currentBackend(r, job.w), &qs, &ss, &ts, r.projectModel.get(), r.riemannModel.get(), r.cubicModel.get(), r.quotientModel.get(), r.tripleModel.get(), r.masterModel.get(), r.qqrhDirectModel.get());
                    }
                    auto text = new std::string(ok ? "HQ saved: " + utf8(job.output.filename().wstring())
                                                   : "HQ cancelled.");
                    PostMessageW(window, NOTICE, 1, reinterpret_cast<LPARAM>(text));
                    continue;
                }
                const std::vector<int> strides =
                    job.singlePass ? std::vector<int>{1}
                    : job.previewOnly ? (options.interaction=="coarse" ? std::vector<int>{16} : std::vector<int>{16,4})
                    : options.refine == "fast"   ? std::vector<int>{16, 4, 1}
                    : options.refine == "direct" ? std::vector<int>{16, 1}
                                                   : std::vector<int>{8, 4, 2, 1};
                for (int stride : strides) {
                    if (generation.load() != job.generation)
                        break;
                    // Publication is atomic; context ownership begins here,
                    // after the compiler has finished all device-only work.
                    accelerated=gpuBackend.renderer();
                    if(accelerated)accelerated->setGpuTuning(options.gpuTuning);
                    cpuPreviewEngine.setThreadLimit(gpuBackend.loading()?compilingThreads:previewThreads);
                    auto f = std::make_unique<Frame>();
                    f->generation = job.generation;
                    f->w = (job.w + stride - 1) / stride;
                    f->h = (job.h + stride - 1) / stride;
                    f->targetW = job.w;
                    f->targetH = job.h;
                    f->stride = stride;
                    f->state = job.state;
                    f->inputAt = job.inputAt;
                    f->projectModel = r.projectModel;
                    f->riemannModel = r.riemannModel;
                    f->cubicModel = r.cubicModel;
                    f->quotientModel = r.quotientModel;
                    f->tripleModel = r.tripleModel;
                    f->masterModel = r.masterModel;
                    f->qqrhDirectModel = r.qqrhDirectModel;
                    const bool cpuPreview=!options.cpu && r.state.family==20 && stride>1 &&
                        (options.previewDevice=="cpu" || (options.previewDevice=="auto" &&
                            std::uint64_t(f->w)*std::uint64_t(f->h)<=65536));
                    Renderer& frameEngine=(!accelerated||rhPreferCpu(options,r.state)||cpuPreview)?cpuPreviewEngine:*accelerated;
                    auto start = Clock::now();
                    if (!frameEngine.render(r, f->w, f->h, double(job.h) / job.w, 1, f->pixels, &generation,
                                       job.generation, stride==1 ? 0 : stride>=8 ? 17 : 19))
                        break;
                    f->ms = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
                    f->gpuMetrics = frameEngine.rhGpuMetrics();
                    f->backend = frameEngine.currentBackend(r, f->w);
                    f->quantumStats = frameEngine.quantumStats();
                    f->sphereStats = frameEngine.sphericalStats();
                    f->spinStats = frameEngine.spinStats();
                    if (!PostMessageW(window, FRAME, 0, reinterpret_cast<LPARAM>(f.get())))
                        break;
                    f.release();
                }
            } catch (const std::exception &e) {
                auto s = new std::string("Render error: " + std::string(e.what()));
                PostMessageW(window, NOTICE, exporting ? 1 : 0, reinterpret_cast<LPARAM>(s));
            }
            if(!exporting)
                PostMessageW(window,RENDER_DONE,static_cast<WPARAM>(job.generation),0);
        }
    }
    void fill(HDC dc, Box b, COLORREF c) {
        HBRUSH brush = CreateSolidBrush(c);
        RECT r = b.rect();
        FillRect(dc, &r, brush);
        DeleteObject(brush);
    }
    void border(HDC dc, Box b, COLORREF c) {
        HBRUSH brush = CreateSolidBrush(c);
        RECT r = b.rect();
        FrameRect(dc, &r, brush);
        DeleteObject(brush);
    }
    void text(HDC dc, Box b, const std::string &s, COLORREF c, UINT align = DT_LEFT) {
        SetTextColor(dc, c);
        SetBkMode(dc, TRANSPARENT);
        RECT r = b.rect();
        auto ws = wide(s);
        DrawTextW(dc, ws.c_str(), int(ws.size()), &r,
                  align | DT_SINGLELINE | DT_VCENTER | DT_END_ELLIPSIS | DT_NOPREFIX);
    }
    void line(HDC dc, int x, int y, int xx, int yy, COLORREF c, int thickness = 1) {
        HPEN p = CreatePen(PS_SOLID, sc(thickness), c);
        auto old = SelectObject(dc, p);
        MoveToEx(dc, x, y, nullptr);
        LineTo(dc, xx, yy);
        SelectObject(dc, old);
        DeleteObject(p);
    }
    void button(HDC dc, Box b, const std::string &s, bool chosen = false, bool on = true) {
        fill(dc, b, chosen ? RGB(69, 49, 104) : RGB(29, 31, 41));
        border(dc, b, chosen ? RGB(184, 156, 244) : RGB(69, 77, 93));
        b.x += sc(4);
        b.w -= sc(8);
        text(dc, b, s,
             !on      ? RGB(103, 110, 129)
             : chosen ? RGB(251, 244, 255)
                      : RGB(217, 221, 234),
             DT_CENTER);
    }
    void slider(HDC dc, Slider s) {
        bool on = enabled(s.id);
        text(dc, {s.box.x + sc(2), s.box.y, s.box.w - sc(4), sc(20)}, label(s.id),
             on ? RGB(224, 229, 241) : RGB(116, 124, 146));
        int left = s.box.x + sc(8), w = s.box.w - sc(16), top = s.box.y + sc(25);
        int knob = left + int(std::clamp((value(s.id) - s.lo) / (s.hi - s.lo), 0., 1.) * w);
        fill(dc, {left, top - sc(2), w, sc(5)}, on ? RGB(62, 69, 81) : RGB(36, 39, 51));
        if (on) {
            fill(dc, {left, top - sc(2), knob - left, sc(5)}, RGB(171, 133, 242));
            fill(dc, {knob - sc(5), top - sc(7), sc(11), sc(15)},
                 active == s.id ? RGB(255, 221, 138) : RGB(229, 225, 247));
        }
    }
    void marker(HDC dc, Box area, Camera cam, C z, COLORREF color) {
        double x = area.x + area.w * (.5 + (z.real() - cam.x) / cam.span),
               y = area.y + area.h * .5 - (z.imag() - cam.y) * area.w / cam.span;
        if (x < area.x + 7 || y < area.y + 7 || x > area.x + area.w - 7 || y > area.y + area.h - 7)
            return;
        line(dc, int(x) - 6, int(y), int(x) + 7, int(y), RGB(4, 4, 8), 3);
        line(dc, int(x), int(y) - 6, int(x), int(y) + 7, RGB(4, 4, 8), 3);
        line(dc, int(x) - 6, int(y), int(x) + 7, int(y), color);
        line(dc, int(x), int(y) - 6, int(x), int(y) + 7, color);
    }
    void drawTrail(HDC dc, Box area, Camera cam, const std::vector<Q> &list, COLORREF col) {
        int saved = SaveDC(dc);
        IntersectClipRect(dc, area.x, area.y, area.x + area.w, area.y + area.h);
        POINT old{};
        bool prev = false;
        for (size_t i = 0; i < list.size(); i++) {
            C z;
            if (!affine(list[i], z)) {
                prev = false;
                continue;
            }
            double x = area.x + area.w * (.5 + (z.real() - cam.x) / cam.span),
                   y = area.y + area.h * .5 - (z.imag() - cam.y) * area.w / cam.span;
            if (!std::isfinite(x) || !std::isfinite(y) || std::abs(x) > 50000 || std::abs(y) > 50000) {
                prev = false;
                continue;
            }
            POINT p{LONG(x), LONG(y)};
            if (prev && std::abs(p.x - old.x) < area.w * 2 && std::abs(p.y - old.y) < area.h * 2)
                line(dc, old.x, old.y, p.x, p.y, col);
            fill(dc, {p.x - sc(2), p.y - sc(2), sc(4), sc(4)}, i == 0 ? RGB(255, 235, 178) : col);
            old = p;
            prev = true;
        }
        RestoreDC(dc, saved);
    }
    void drawLegacyLambdaGuides(HDC dc) {
        if(!legacyLambdaMode()) return;
        const auto b=legacyLambdaBox();
        const auto &cam=state.cameras[2];
        const int saved=SaveDC(dc);
        IntersectClipRect(dc,b.x,b.y,b.x+b.w,b.y+b.h);
        const double axisX=b.x+b.w*(.5-cam.x/cam.span),axisY=b.y+b.h*(.5+cam.y/cam.span);
        if(std::isfinite(axisX) && axisX>=b.x && axisX<b.x+b.w)
            line(dc,int(axisX),b.y,int(axisX),b.y+b.h,RGB(111,112,135));
        if(std::isfinite(axisY) && axisY>=b.y && axisY<b.y+b.h)
            line(dc,b.x,int(axisY),b.x+b.w,int(axisY),RGB(111,112,135));
        bool previous=false;int oldX=0,oldY=0;
        for(int i=0;i<=192;++i) {
            const double angle=2*pi*i/192;
            const double x=b.x+b.w*(.5+(std::cos(angle)-cam.x)/cam.span);
            const double y=b.y+b.h*(.5-(std::sin(angle)-cam.y)/cam.span);
            if(!std::isfinite(x) || !std::isfinite(y) || std::abs(x)>50000 || std::abs(y)>50000) {previous=false;continue;}
            if(previous)line(dc,oldX,oldY,int(x),int(y),RGB(198,177,220));
            oldX=int(x);oldY=int(y);previous=true;
        }
        marker(dc,b,cam,C(1,0),RGB(255,221,138));
        text(dc,{b.x+sc(9),b.y+sc(5),b.w-sc(18),sc(20)},connectedLambdaMode() ? "Im outer lambda / unit circle guide" : "Im lambda / unit circle guide",RGB(233,228,241));
        text(dc,{b.x+sc(9),b.y+b.h-sc(25),b.w-sc(18),sc(20)},connectedLambdaMode() ? "Re outer lambda" : "Re lambda",RGB(233,228,241),DT_RIGHT);
        RestoreDC(dc,saved);
    }
    void drawOrbit(HDC dc) {
        if (!orbitVisible || state.family == 9 || sampleMode() || legacyLambdaMode())
            return;
        if(directAtlasMode()) {
            if(state.plane==1) drawTrail(dc,{0,0,width,canvasH},state.cameras[1],trailA,RGB(255,222,148));
            return;
        }
        if (state.plane == 0) {
            drawTrail(dc, {0, 0, width, canvasH}, state.cameras[0], trailA, RGB(255, 222, 148));
            return;
        }
        marker(dc, {0, 0, width, canvasH}, state.cameras[state.plane],
               state.plane == 1 ? state.c : state.alpha, RGB(255, 221, 138));
        Box b{width - sc(296), sc(16), sc(280), std::min(sc(190), canvasH - sc(30))};
        if (b.h < sc(95))
            return;
        fill(dc, b, RGB(14, 15, 25));
        border(dc, b, RGB(146, 111, 204));
        text(dc, {b.x + sc(10), b.y + sc(5), b.w - sc(20), sc(20)}, "Dynamical orbit / z plane",
             RGB(236, 222, 255));
        Box plot{b.x + sc(10), b.y + sc(28), b.w - sc(20), b.h - sc(50)};
        Camera cam = (state.family == 6 || state.family == 7 || state.family == 8 || state.family == 10) ? state.cameras[0] : state.family ? Camera{state.family == 2 ? .5 : 0, 0, 5.4} : Camera{2.1, 0, 9};
        line(dc, plot.x, plot.y + plot.h / 2, plot.x + plot.w, plot.y + plot.h / 2, RGB(42, 44, 59));
        if (!state.family) {
            marker(dc, plot, cam, {2, 0}, RGB(145, 148, 171));
            marker(dc, plot, cam, {4, 0}, RGB(145, 148, 171));
        }
        drawTrail(dc, plot, cam, trailA, RGB(255, 221, 138));
        drawTrail(dc, plot, cam, trailB, RGB(171, 151, 247));
        text(dc, {b.x + sc(10), b.y + b.h - sc(22), b.w - sc(20), sc(19)},
             state.family ? "gold: selected seed orbit" : "gold: +critical / lilac: -critical",
             RGB(153, 160, 183));
    }
    void hilbertInspector(HDC dc) {
        if(!hilbertField() || !hilbertInspectorVisible || help || width<sc(520) || canvasH<sc(270))return;
        const auto& model=qqrh_hilbert::cached(state.alpha.real());
        const double overlap=model.coherence(hilbertPicked);
        const double distance=model.projective_distance(hilbertPicked);
        const double metric=model.projective_metric(hilbertPicked);
        marker(dc,{0,0,width,canvasH},state.cameras[state.plane],hilbertPicked,RGB(255,222,148));
        Box b{width-sc(330),sc(14),sc(316),sc(223)};
        fill(dc,b,RGB(12,18,30));
        border(dc,b,RGB(136,155,187));
        auto old=SelectObject(dc,smallFont);
        text(dc,{b.x+sc(12),b.y+sc(7),b.w-sc(24),sc(19)},
            "HILBERT INSPECTOR / I: hide",RGB(233,219,168));
        text(dc,{b.x+sc(12),b.y+sc(30),b.w-sc(24),sc(18)},
            "z = "+complexText(hilbertPicked),RGB(210,220,237));
        text(dc,{b.x+sc(12),b.y+sc(51),b.w-sc(24),sc(18)},
            "Squared overlap = "+num(overlap,7),RGB(210,220,237));
        text(dc,{b.x+sc(12),b.y+sc(72),b.w-sc(24),sc(18)},
            "Projective angle = "+num(distance,6)+" rad",RGB(210,220,237));
        text(dc,{b.x+sc(12),b.y+sc(93),b.w-sc(24),sc(18)},
            "Local metric g = "+num(metric,7),RGB(210,220,237));
        // Log weights avoid overflow when inspecting far from the origin.
        std::array<double,qqrh_hilbert::degree+1> weights{};
        const double scale=std::max(std::abs(hilbertPicked.real()),std::abs(hilbertPicked.imag()));
        if(scale>0. && std::isfinite(scale)) {
            const C reduced=hilbertPicked/scale;
            const double logMagnitude=std::log(scale)+.5*std::log(std::norm(reduced));
            double largest=-std::numeric_limits<double>::infinity();
            for(int n=0;n<=qqrh_hilbert::degree;++n) {
                weights[n]=2.*std::log(std::abs(model.v[n]))+2.*n*logMagnitude;
                largest=std::max(largest,weights[n]);
            }
            double total=0.;
            for(double& value:weights){value=std::exp(value-largest);total+=value;}
            for(double& value:weights)value/=total;
        } else weights[0]=1.;
        const int plotX=b.x+sc(13),plotY=b.y+sc(146),plotW=b.w-sc(26),plotH=sc(44);
        text(dc,{plotX,b.y+sc(119),plotW,sc(18)},
            "Coefficient probabilities |v_n z^n|^2 / K",RGB(155,177,207));
        line(dc,plotX,plotY+plotH,plotX+plotW,plotY+plotH,RGB(71,87,110));
        for(int n=0;n<=qqrh_hilbert::degree;++n) {
            const int left=plotX+n*plotW/(qqrh_hilbert::degree+1);
            const int right=plotX+(n+1)*plotW/(qqrh_hilbert::degree+1);
            const int barH=int(std::round(plotH*weights[n]));
            if(barH>0)fill(dc,{left+sc(2),plotY+plotH-barH,std::max(1,right-left-sc(4)),barH},RGB(210,174,94));
            text(dc,{left,b.y+sc(195),right-left,sc(17)},std::to_string(n),RGB(154,173,199),DT_CENTER);
        }
        SelectObject(dc,old);
    }
    void helpOverlay(HDC dc) {
        if (!help)
            return;
        Box b{sc(18), sc(18), std::min(sc(775), width - sc(36)), std::min(sc(334), canvasH - sc(32))};
        if (b.h < sc(190))
            return;
        fill(dc, b, RGB(16, 17, 27));
        border(dc, b, RGB(184, 156, 244));
        std::vector<std::string> lines = {
            "FRAKTALFAKTUR / native Windows controls",
            "Drag: pan   Wheel: cursor-centred zoom   Right-click: pick + orbit",
            "J / C / A: z-start, c, alpha plane   1-4: presets   F11: fullscreen",
            "E: explore picked parameter   B: back   R: reset camera",
            "S: save BMP   Ctrl+S: save view   Ctrl+O: load view   H: help",
            "Esc: cancel HQ, dismiss help, or leave fullscreen",
            "Click sliders, drag, or use the wheel for fine adjustments.",
            "R(z) = c + (3/5) alpha z(8-3z) / ((2-z)(4-z))",
            "Critical seeds: (12 +/- 4i)/5. Poles: 2 and 4.",
            "Poles go to infinity; infinity returns to c - (9/5) alpha.",
            "Period labels are numerical detections up to 16, not proofs.",
            "GPU: float32. Deep zoom switches to double-precision CPU.",
            "BMP saves the finished canvas. HQ uses 2 x 2 samples per pixel."};
        if (state.family)
            lines = {"FRAKTALFAKTUR / " + familyName(state.family) + " generation " +
                         std::to_string(state.generation),
                     "Drag: pan   Wheel: zoom   Right-click: pick + orbit   Tab: "
                     "catalogue",
                     "J/C/A: z-start / first parameter / second parameter   1-4: presets",
                     "Newton: z <- z - f/f'     Halley: z <- z - 2ff'/(2f'^2 - ff'')",
                     state.family == 1 ? "P[0]=1; P[n+1]=1+a*z*P[n](z^2)+b*z^4*P[n](z^4)"
                                       : "Graph: " + std::to_string(state.generation) +
                                             " x 1 x 1 cells, source (0,0,0) to (L,1,1)",
                     state.family == 1 ? "f(z)=(1-z)P[n](z); a and b are complex branch weights."
                                       : "f(p)=C[L](p)-target; complex relaxation multiplies the step.",
                     state.family == 1 ? "Generation is selected in the catalogue."
                                       : "C[L] counts connected bond subsets exactly; "
                                         "rendering is numerical.",
                     "Root phase: hue by final root, bands by approach time.",
                     "Dark purple: iteration limit; amber: derivative pole; grey: cutoff.",
                     "Convergence uses a scaled relative polynomial residual.",
                     "S: BMP   Ctrl+S: view   Ctrl+O: load view   F11: fullscreen",
                     "Finite diagnostics; no mathematical convergence proof is claimed.",
                     "Deep zoom uses double-precision CPU. HQ uses 2 x 2 samples."};
        if (state.family == 3) {
            lines[4] = "S(n)=(3n+1)/2^a; words (a1,...,aN) use a in {1,2,3}.";
            lines[5] = "Each word fixes n=r mod 2^(A+1), where A=sum(a).";
            lines[6] = "Embed: rho(q)=sum[j=1..A] (2*bit_j(r)-1)*q^(j-1).";
            lines[7] = "f(z)=product(z-rho); complex relaxation multiplies the step.";
            lines[8] = "Cylinder basins: colour by reached word; dark = iteration limit.";
            lines[9] = "Convergence: nearest-root distance / (1+root magnitude) < 3e-6.";
            lines[11] = "The complex embedding is a chosen geometry for exact "
                        "residue cylinders.";
        }
        if (state.family == 4) {
            auto &t = diagnostic.gardenData->transfer;
            lines[4] = "Odd states 1.." + std::to_string(t.cutoff) + "; follow S(n)=(3n+1)/2^valuation.";
            lines[5] = "At 1 or an outgoing boundary jump: restart uniformly, with feedback q.";
            lines[6] = "T(q)=Q+q*u*b^T; det(zI-T)=z^(M-L)*D(z).";
            lines[7] = "D(z)=z^L-q*sum[l=1..L] (h_l/M)*z^(L-l).";
            lines[8] = "h_l counts states taking exactly l steps to the next restart.";
            lines[9] = "This garden iterates D; the fixed zero factor is removed.";
            lines[11] = "Finite transfer model; q=1 is stochastic. Complex q deforms its spectrum.";
        }
        if (state.family == 5) {
            lines[4] = "Cayley core: q(z)=4z/(1+z)^2; q=1-((z-1)/(z+1))^2.";
            lines[5] = "P0=z; Q0=1; Pnext=4PQ; Qnext=(P+Q)^2.";
            lines[6] = "f=P[depth]-target*Q[depth]; maximum degree 2^depth.";
            lines[7] = "Newton/Halley acts on this polynomial; relaxation multiplies its step.";
            lines[9] = "Common jet scaling avoids overflowing powers; poles are not roots.";
            lines[11] = "Chosen complex dynamics from the exact kernel identity.";
        }
        if (state.family == 6) {
            lines[0] = "FRAKTALFAKTUR / " + lfunctionName(state.generation);
            lines[4] = "f(s)=L(s,chi)-target; target=0 selects the L-function zeros.";
            lines[5] = "Euler-Maclaurin continuation and functional equation; analytic jets.";
            lines[6] = "Re s in [-8,12], Im s in [-80,80]; grey means evaluator cutoff.";
            lines[7] = "CPU double precision. Converged basins use absolute L residual.";
            lines[9] = "Complex target changes the level set; relaxation changes the iteration.";
            lines[11] = "A numerical exploration of the RH/GRH objects; not a zero certificate.";
            lines[12] = "See RH_GRH_MATHEMATICS.txt and LFUNCTION_NUMERICS.txt.";
        }
        if (state.family == 7) {
            lines[0] = "FRAKTALFAKTUR / Primitive GRH Cayley / chi mod 5";
            lines[4] = "w=tanh(log(2)*(s-1/2)); sheet inverse adds i*pi*k before division by log(2).";
            lines[5] = "At target=0: f=A(s)L(s,chi), A=(5/pi)^((s+1)/2) Gamma((s+1)/2).";
            lines[6] = "A common completion factor is divided from f,f',f''; Newton/Halley is unchanged.";
            lines[7] = "GRH carrier: Re(s)=1/2 maps exactly to Re(w)=0, drawn as the vertical guide.";
            lines[8] = "Sheet k=" + std::to_string(state.generation) + " contains the selected vertical branch of atanh.";
            lines[9] = "Convergence uses |L-target|, avoiding false zeros from Gamma decay.";
            lines[11] = "This is a numerical phase portrait of the primitive GRH object, not a proof certificate.";
            lines[12] = "See PRIMITIVE_GRH_CAYLEY_MATHEMATICS.txt.";
        }
        if (state.family == 8) {
            lines[0] = "FRAKTALFAKTUR / Primitive Character Pair / chi5 x chibar5";
            lines[4] = "P(s)=Lambda(s,chi5)*Lambda(s,chibar5), with both primitive quartic characters.";
            lines[5] = "At target=0 the attractors are the union of the chi5 and chibar5 zero sets.";
            lines[6] = "The common completion factor A(s)^2 is divided from all jets; Newton/Halley is unchanged.";
            lines[7] = "L(s,chibar)=conjugate(L(conjugate(s),chi)); analytic derivative jets use the same identity.";
            lines[8] = "The paired basins expose the conjugation / functional-equation mirror geometry directly.";
            lines[9] = "A nonzero target is a normalized exploratory deformation of the pair product.";
            lines[11] = "Numerical phase portrait only; convergence is not a GRH certificate.";
            lines[12] = "See PRIMITIVE_PAIR_MATHEMATICS.txt.";
        }
        if (state.family == 9) {
            lines = {"FRAKTALFAKTUR / literal GRH Cayley observable diagnostic",
                     "Shared plane: x=log2(a), y=logit(r); drag/zoom moves all three panels together.",
                     "LEFT: phase of Theta_chi(z+), brightness from log|Theta|.",
                     "MIDDLE: phase + magnitude of R_chi(z+)=Theta'/Theta.",
                     "RIGHT: signed magnitude of C_chi(a,r), the literal two-point observable.",
                     "z+ = a-2ar+2ia sqrt(r(1-r)); z- = conjugate(z+).",
                     "Theta_chi(z)=completedPair_chi(1/2+sqrt(z)) on the principal slit branch.",
                     "C=4a[z-R(z-)-z+R(z+)]/(z--z+).",
                     "Critical-line zeros approach the slit boundary r -> 1; off-line zeros would be interior failures.",
                     "Right-click reports numerical a,r and C at the clicked point.",
                     "1-4: boundary/wide/high/extreme-boundary presets. Palette buttons change contrast only.",
                     "CPU double precision; all three panels share each expensive completed-L evaluation.",
                     "Numerical visualization only. See GRH_OBSERVABLE_MATHEMATICS.txt."};
        }
        if (state.family == 10) {
            lines = {"FRAKTALFAKTUR / GRH Phase-Feedback Garden / primitive chi mod 5",
                     "Cayley plane w on sheet 3; the completed primitive-L zeros remain fixed attractors.",
                     "Base motion is exact completed-L Newton or Halley after common-factor cancellation.",
                     "A bounded residual-phase kick is added: kappa*Phi(f;tau)/(1+0.4|w|^2).",
                     "Phi=tanh(0.8|f|)*(f/|f|)*exp(i*tau*arg(f)); Phi tends to zero at every zero.",
                     "Therefore feedback changes basin boundaries without moving the GRH zero attractors.",
                     "kappa controls turbulence; tau twists residual phase; complex relaxation rotates the root step.",
                     "Preset 1: flowing default. 2: eddy zoom. 3: turbulent weave. 4: Halley feedback.",
                     "The vertical Cayley carrier Re(w)=0 is the image of Re(s)=1/2.",
                     "This nonlinear feedback map is exploratory dynamics built from the formal GRH object.",
                     "It is NOT itself a new GRH equivalence or proof certificate.",
                     "CPU double precision; the launcher accepts --threads and --refine options.",
                     "See GRH_PHASE_FEEDBACK_MATHEMATICS.txt."};
        }
        if (state.family == 11) {
            lines = {"FRAKTALFAKTUR / Spectral Jet Determinant",
                     "F(x)=x^2+c+eta*d^2*(x^2-d^2)/(x^2+d^2)^2; direct rational iteration.",
                     "Two equal positive atoms: q=1/2-d and 1/2+d; centered coordinate x=z-1/2.",
                     "J: x-start plane. C: complex c plane. A: complex eta plane.",
                     "The parameter planes follow the selected seed, initially x=0.",
                     "Right-click in J traces an orbit and selects its start as the parameter seed.",
                     "Pole markers toggles the guides at +id and -id. Seed=0 restores the marked seed.",
                     "Drag sliders or use the wheel. Spectral spacing d controls the two-atom separation.",
                     "1: lace. 2: islands. 3: filaments. 4: quadratic control (eta=0).",
                     "Escape palette: dark means not escaped at this budget; magenta means numerical failure.",
                     "Pole trap and orbit phase colour finite orbit geometry, including unresolved points.",
                     "S: BMP; Ctrl+S/O: save/load view; HQ: 2x2 samples; F11: fullscreen; Tab: catalogue.",
                     "GPU at normal zoom; double-precision CPU at deep zoom. See SPECTRAL_JET_MATHEMATICS.txt."};
        }
        if (spectralGen2()) {
            lines = {"FRAKTALFAKTUR / Spectral Jet Generation 2: three atoms",
                     "F(x)=x^2+c+eta*(R_s*R_(s+2)/R_(s+1)^2-1); choose jet order s=1,2,3.",
                     "Centered atoms: -d, 0.35d, +d. Their positive locations are q=1/2+atom.",
                     "Weights: (1-m)/2, m, (1-m)/2. Middle weight m changes the spectral geometry.",
                     "J: initial x. C: complex c. A: complex eta. Parameter planes follow the selected seed.",
                     "Right-click in J traces an orbit and selects the complex parameter seed.",
                     "Atom guides mark the three real atoms, not the map's complex poles.",
                     "Sliders: c, eta, spacing d, middle weight m, jet order s and both seed components.",
                     "1: Triad. 2: Split branches. 3: Second jet. 4: Third jet.",
                     "Smooth escape: dark is unresolved at the budget; magenta is numerical failure.",
                     "Denominator trap measures cancellation near R_(s+1)=0, not Euclidean pole distance.",
                     "S: BMP; Ctrl+S/O: views; HQ: 2x2 samples; F11: fullscreen; Tab: generations.",
                     "Canceled polynomial evaluation crosses spectral atoms safely. See SPECTRAL_GENERATION_2.txt."};
        }
        if (quantumMode()) {
            lines={"FRAKTALFAKTUR / Quantum Measurement Histories",
                "Each point is a conditioned weak-measurement state on the Bloch sphere.",
                "K_i=(I+epsilon*n_i.sigma)/sqrt(N*(1+epsilon^2)); sum K_i^dag K_i=I.",
                "Rotate the state by the unitary kick, choose a Born-weighted outcome, normalize.",
                "J: rotatable sphere. C: equal-area longitude/height atlas. A: mean/history comparison.",
                "Drag: rotate. Wheel: zoom. Arrow keys: rotate. R: reset orientation and zoom.",
                "1: tetrahedron. 2: octahedron. 3: icosahedron. 4: kicked dodecahedron.",
                "Strength changes the measurement. Geometry changes its directions. Kick changes precession.",
                "New history / N changes the replay seed. Reusing a seed reproduces the same chains.",
                "Mean view: cyan is exact ensemble decay, amber is a finite-chain estimate (before burn-in).",
                "Density uses retained states after burn-in; finite samples are not a dimension certificate.",
                "S: BMP; Ctrl+S/O: views; HQ: up to 4x samples; Tab: catalogue; F11: fullscreen.",
                "Parallel CPU sampling; cached points support camera and colour changes without resampling."};
        }
        if (sphereMode()) {
            lines={"FRAKTALFAKTUR / Spherical Fractals / "+std::to_string(state.generation),
                state.generation==3 ? "McMullen: f(z)=z^3+lambda/z^3 on the Riemann sphere."
                                    : "Four circle reflections generate a spherical group limit set.",
                state.generation==3 ? "Backward iteration chooses uniformly among the six preimages, with multiplicity."
                                    : "Reduced reflection words choose among the other three generators at each step.",
                state.generation==3 ? "Lambda magnitude and phase change the map; topology depends on the parameter."
                                    : "Tangency gives an Apollonian gasket; separating the circles gives Schottky dust.",
                "J: rotatable sphere. C: equal-area longitude/height atlas. A: front/back hemispheres.",
                "Drag: rotate. Wheel: zoom. Arrow keys: rotate. R: reset orientation and zoom.",
                "1-4: presets. Circle separation applies to Schottky; lambda controls apply to McMullen.",
                "Conformal display boost redistributes the visible geometry without resampling the dynamics.",
                "N: new reproducible sample. E/B: raise/lower the sample budget. Replay seed saves with views.",
                "Density measures sampled words or inverse branches, not Born probabilities or uniform area.",
                "Finite samples reveal geometry; they do not certify dimension or parameter-specific topology.",
                "S: BMP; Ctrl+S/O: views; HQ: up to 4x samples; Tab: catalogue; F11: fullscreen.",
                "Parallel CPU sampling; cached points support camera and colour changes without resampling."};
        }
        if (spinMode()) {
            lines={"FRAKTALFAKTUR / Dissipative Spin Chaos / "+std::to_string(state.generation),
                "Each step: tilt about y, twist about z by k*height, then damp toward the south pole.",
                "The spin length stays one. Initial points are seeded; subsequent dynamics is deterministic.",
                state.generation==3 ? "J: parameter atlas. C: selected sphere. A: selected equal-area atlas."
                                    : "J: sphere. C: equal-area atlas. A: front/back hemispheres.",
                state.generation==3 ? "Atlas: drag pans, wheel zooms at cursor, right-click selects k/tau, E opens its sphere."
                                    : "Drag: rotate. Wheel: zoom. Arrows: rotate. R: reset camera. E/B: sample budget.",
                state.generation==3 ? "B returns to the parameter atlas. R resets its bounds. Grid width controls sampling resolution."
                                    : "1-4: presets. N selects new initial points. Camera and palette changes reuse cached trajectories.",
                "Cloud colours: visit density, finite-time tangent growth, or mean log area change.",
                "Atlas colours: finite-time tangent growth, mean log area change, or final spin height.",
                "Burn-in discards transients; cloud burn-in is at least the Lyapunov averaging window.",
                "The stretching estimate propagates an analytic tangent derivative along each orbit.",
                "Finite samples and positive finite-time growth do not certify chaos, stationarity or dimension.",
                "S: BMP; Ctrl+S/O: views; HQ: denser cloud or parameter grid; Tab: catalogue; F11: fullscreen.",
                "Parameter scans use the same seeded initial state in each cell; coexisting basins can be missed."};
        }
        if(projectMode()) {
            lines={projectTitle(state.generation),
                "Drag: pan   Wheel: cursor-centred zoom   Right-click: trace orbit   1-4: presets",
                "This entry uses the complex initial-state plane. Its update rule is fixed.",
                formula(state),
                state.generation==1 ? "GRH: finite polynomial model of the centered completed chi5/chibar5 pair."
                : state.generation==2 ? "PI: actual binary pairs at n=1024+offset+t, minus the fixed-anchor drift 1024*c(n)."
                : state.generation==3 ? "NSE: spectral energies of an exact periodic shear; viscosity changes coefficients."
                : state.generation==4 ? "COLLATZ: the entire shortcut extension is applied twice per displayed step."
                : "P_c(d,q): q=2 GP-4 quintic in the complex dimension variable sigma=2d-1.",
                state.generation==1 || state.generation==4 ? "Direct dynamics: finite cycle detection and unresolved orbit diagnostics."
                                                         : "Newton dynamics: root basins with numerical convergence thresholds.",
                "An iteration limit or large-magnitude cutoff is not an escape theorem.",
                "These complex dynamical maps do not prove the underlying project conjectures.",
                "Palettes: basin/cycle phase, convergence/recurrence time, finite orbit structure.",
                "S: BMP   Ctrl+S: view   Ctrl+O: load view   R: reset camera   Tab: catalogue",
                "The command --plot-five DIRECTORY saves all five maps and a labelled contact sheet.",
                "CPU double precision. HQ uses 2 x 2 samples per pixel."};
        }
        if(riemannMode()) {
            lines={riemannTitle(state.generation),
                "Drag: pan   Wheel: cursor-centred zoom   Right-click: trace orbit   1-4: presets",
                state.generation>=4 ? "Five fixed moments c0..c4 determine this polynomial; the moment order is H2."
                                    : "Finite witness H1: two nodes. H2: three nodes with one prescribed anchor.",
                formula(state),
                state.generation==7 ? "Joukowski coordinate: t=(z+1/z)/2, with t=(q-mu)/sigma; the cleared expression is polynomial."
                    : state.generation>=4 ? "Standardized coordinate: z=(q-mu)/sigma, with mean mu and positive scale sigma."
                                          : "kappa(z)=(27/4) z (1-z)^2; generations differ in their polynomial pullback.",
                "Newton: z <- z - alpha P/P'. Halley: z <- z - alpha 2PP'/(2P'^2-PP'').",
                "The complex relaxation alpha is nonzero and has magnitude at most 2.",
                state.generation==6 ? "Jacobi closure sets the terminal diagonal to zero; this is an explicit finite model choice."
                                    : "The witness uses finite rounded RH-branch moments; its construction is numerical.",
                state.generation>=4 ? "The witness anchor is diagnostic only; these fixed moment polynomials do not depend on it."
                                    : "H2 anchor changes may fail positivity checks; model status reports such failures.",
                "Basin colours and recurrence diagnostics do not prove the infinite RH binding.",
                state.generation>=4 ? "BMP exports save moments, standardization, witness diagnostics and the actual iteration."
                                    : "BMP exports save the actual model's moments, nodes, weights and residual.",
                "S: BMP   Ctrl+S: view   Ctrl+O: load view   R: reset camera   Tab: catalogue",
                "CPU double precision. HQ uses 2 x 2 samples per pixel."};
        }
        if(cubicMode()) {
            lines={"RH / Jacobi cubic Julia atlas",
                "Direct iteration: f_c(z)=z^3-A*z+c, A=1+h+b^2/3 from five fixed RH moments.",
                "The finite Jacobi characteristic cubic is translated to its depressed cubic form.",
                "Atlas: each pixel is c; follow BOTH critical points +/- sqrt(A/3).",
                "Right-click the atlas to select c and open its Julia plane. C or B returns to the atlas.",
                "J: selected Julia plane. Right-click there to trace an initial-state orbit.",
                "An escape-radius crossing detects escape numerically; budget-limited orbits stay unresolved.",
                "Nonescape within a finite budget does not certify a connected Julia set.",
                "1: critical atlas. 2: symmetric Julia. 3: critical-fixed Julia. 4: complex Julia.",
                "Palettes: critical fates, escape time, finite orbit structure. A is fixed; c changes the map.",
                "Drag: pan. Wheel: zoom. R: reset camera. Tab: catalogue. H: help.",
                "S: BMP; Ctrl+S/O: views; HQ: 2 x 2 samples per pixel. CPU double precision.",
                "The finite RH moment model and this cubic dynamics do not prove RH."};
        }
        if(quotientMode()) {
            lines={"RH / Jacobi quotient Julia atlas",
                "Direct rational iteration: f_c(z)=lambda*(z-b-1/z)+c, lambda=1/sqrt(h)>1.",
                "Consecutive orthogonal polynomials give (z*z-b*z-1)/(sqrt(h)*z).",
                "Five fixed RH moments determine b and h; this is a finite Jacobi quotient.",
                "Atlas: each pixel is c; follow BOTH critical points +i and -i.",
                "Right-click the atlas to open its Julia plane. C or B returns; J selects Julia.",
                "Julia view colours the infinity basin. Right-click to trace a finite orbit prefix.",
                "Pole 0 maps to infinity; infinity is fixed, with multiplier sqrt(h)<1.",
                "An escape-radius crossing detects the infinity basin numerically; finite nonescape is unresolved.",
                "1: critical atlas. 2: symmetric Julia. 3: critical-fixed Julia. 4: complex Julia.",
                "Palettes: critical fates / escape phase, escape time, finite orbit structure.",
                "Drag: pan. Wheel: zoom. R: reset. Tab: catalogue. S: BMP; Ctrl+S/O: views.",
                "The finite RH moment model and this rational dynamics do not prove RH."};
        }
        if(tripleMode()) {
            lines={momenttriple::title(state.generation),momenttriple::formula(state.generation),
                "A=1+h+b*b/3; a=1+b*b/4; eta=sqrt(h), from the fixed RH point moments.",
                state.generation==1 ? "Conjugation reverses orientation; the second iterate is holomorphic of degree nine." :
                state.generation==2 ? "Q+(z)=z*z-a+c, then Q-(z)=z*z-a-c. A displayed step is the full return." :
                                     "The rational map has degree three; its simple pole 0 maps to infinity.",
                "Atlas: each pixel is c; every finite critical orbit is followed.",
                "Right-click atlas: choose c and open Julia. J: Julia. C / B: atlas.",
                "Right-click Julia: trace a finite orbit. Dark nonescape stays unresolved.",
                "Palette 0 encodes the escaped-critical subset; it is not a connectedness label.",
                "1: atlas. 2: c=0 Julia. 3: critical fixed. 4: complex Julia.",
                "Drag: pan. Wheel: zoom. R: reset. Tab: catalogue. H: help.",
                "S: BMP; Ctrl+S/O: views. HQ: higher resolution. CPU double precision.",
                "These three iterated maps are visualization constructions from finite moments, not RH proofs."};
        }
        if(extraMode()) {
            lines={std::string(frk_extra::extra_catalogue[state.generation-1].label),formula(state),
                "Fractal dynamics added to the original 18-map RH Lab GPU v1.8 catalogue.",
                "J: starting-z plane. C/B: parameter-c portrait. P: cycle palette. 1-4: presets.",
                "Parameter portraits start from z=0; Newton and Halley maps always start from each pixel.",
                "Right-click a parameter to open its Julia plane; right-click Julia to trace finite orbits.",
                "Adjust real lambda and the iteration budget. The two selected-c sliders apply to Julia views.",
                "GPU computes Shader Model 5 fractal pixels; deep zoom automatically uses CPU double precision.",
                "GPU and CPU floating-point classifications may diverge near fractal boundaries and poles.",
                "Dark pixels mean finite-budget nonescape or nonconvergence, not a certified set membership.",
                "Drag pan / wheel zoom / R reset / Tab catalogue / S BMP / Ctrl+S/O views / H help.",
                "This is a mathematical visualization, not a proof of RH, Collatz, or an L-function hypothesis."};
        }
        if(hilbertMode()) {
            lines={std::string(frk_extra::extra_catalogue[state.generation-1].label),formula(state),
                "Finite Hardy space: 13 coefficients, n=0..12; b_n=200^n*a_n/a_0 from the original Theta table.",
                "R is the scaled radius. Physical Theta radius r=200R. v_n=b_n*R^n, normalized to norm 1.",
                "Overlap: |<w(z),v>|^2/||w(z)||^2. Projective ray angle = acos(sqrt(overlap)), in radians.",
                "Metric: g=partial_z partial_conj(z) log ||w(z)||^2; palette displays g/(1+g).",
                "Right-click either field: inspect angle, local metric, and 13 probabilities. I: show/hide inspector.",
                "Kernel H(z)=sum v_n^2 z^n. H=0 means w(z) is orthogonal to v in this finite model.",
                "Kernel Newton follows z-H/H'; kernel feedback follows z^2+c-.35*H/H'.",
                "Shift and angular entries use c_k=<B^k v,v>; their existing dynamics are preserved.",
                "J: Julia / starting-z plane. C/B: parameter portrait for shift or kernel-feedback fractals.",
                "Drag pan / wheel zoom / R reset / Tab catalogue / P palette / S BMP / Ctrl+S/O views.",
                "Finite numerical geometry and new dynamical maps. No full-Theta tail bound or RH proof."};
        }
        if(masterMode()) {
            if(labControls) {
                lines={masterTitle(state.generation),
                    "T: preset bank (Original Lab / Currents / Weaves). 1..4: bank presets. V: main controls.",
                    "F: map. D: Julia chart. I: coloring. K: Taylor degree. P: three original palettes.",
                    "Cayley currents multiply the host by (1+d)/(1-d); conjugate currents reverse orientation.",
                    "Alternating tides switch ordinary / conjugate carriers; inertial eddies add previous-step motion.",
                    "Twist rotates feedback in maps 9 / 11. Blend mixes z with conjugate(z) in maps 10 / 11.",
                    "Charts: Cartesian z=w; reciprocal z=1/w; log-polar z=exp(w); Joukowski z=(w+1/w)/2.",
                    "Charts change Julia pixels only. Orbit traces follow inverse branches; chart poles end a trace.",
                    "Colors: original; endpoint phase; circle trap; rotated cross trap; orbit winding;",
                    "feedback phase; finite-tail period; orbit radius. Exposure and frequency change coloring.",
                    "Power sets the polynomial host; gain scales feedback; relaxation mixes old and proposed states.",
                    "Memory controls the inertial map. Trap radius / angle shape the circle and rotated cross.",
                    "J: Julia. C/B: c. A: lambda. Tab: catalogue. S: BMP; Ctrl+S/O: views. HQ: supersampling.",
                    "GPU float32; --cpu uses double/long-double. Numerical limits stay unresolved."};
            } else {
                lines={masterTitle(state.generation),formula(state),
                    "V: lab controls. T: preset bank. F: map. D: Julia chart. I: coloring. K: Taylor N=8 / 16 / 32.",
                    "J: Julia. C/B: c portrait. A: complex lambda. Portraits follow fixed physical seed z0=0.",
                    "R is the normalized angular radius; level selects the candidate's equipotential or count level.",
                    "Julia charts affect pixels and orbit traces; c and lambda coordinates remain unchanged.",
                    "Right-click a parameter plane to select c or lambda and open Julia; Julia click traces an orbit.",
                    "Main presets 1..4: default Julia; seeded c portrait; degree 8 Julia; degree 32 Julia.",
                    "Lab banks: Original Lab / Currents / Weaves. V exposes their presets and Twist / Blend sliders.",
                    "P: three original palettes. [ / ]: hue. Comma / period: bands. Drag: pan. Wheel: zoom.",
                    "R: reset camera. Tab: catalogue. H: help. F11: fullscreen. Save orbit: displayed chart trace.",
                    "S: BMP; Ctrl+S/O: views. HQ: 2 x 2 samples. Exports include all experimental parameters.",
                    "Magnitude |z|>16 is a display cutoff; finite-tail periods are numerical recurrence observations.",
                    "GPU float32; --cpu: CPU reference. Finite budgets and numerical limits do not certify RH."};
            }
        }
        if(legacyLambdaMode()) {
            lines={"Original BKQR lambda map / native catalogue",formula(state),
                "This entry scans complex lambda using the original finite ensemble and host-map policies.",
                "A: lambda plane only. F: host. B: GREEDY / SUB / ADD. C: RAW / SELFDUAL / COMPACT.",
                "Hosts: Julia, Newton, Tricorn, Burning Ship, Nova. The power and complex c define the host.",
                "P / M: cycle nine metrics. Palette row also selects previous or next metric.",
                "Candidate: gradient of sampled statistics. Lyap: finite separation growth, not a rigorous exponent.",
                "Quality and precision presets choose the finite grid and seed budget; step budget truncates orbits.",
                "The complex domain is a square with equal axes; default bounds are [-2,2] x [-2,2].",
                "Drag inside the square: pan. Wheel: cursor-centred zoom. Right-click: lambda coordinate readout.",
                "1: raw greedy. 2: raw subtract. 3: selfdual greedy. 4: compact greedy.",
                "R: reset bounds. S: BMP. Ctrl+S/O: views. Tab: catalogue. H: help. F11: fullscreen.",
                "These sampled diagnostics reproduce the old dynamics; numerical textures do not certify phase transitions."};
        }
        if(connectedRhLambdaMode()) {
            lines={"BKQR outer lambda / "+std::string(masterTitle(state.generation)),
                "RH host Q_j(z)=z^2+c+eta*response_j(z); outer update z_next=z+(1+sigma*H(lambda))*(Q_j(z)-z).",
                "Lambda is the complex pixel coordinate. Inner complex eta is an independent RH-host control.",
                "B: branch GREEDY / SUB / ADD. C: H = RAW / SELFDUAL / COMPACT. K: Taylor N = 8 / 16 / 32.",
                "Choose Master, Sub_1 or Sub_2 in the catalogue. Sub_2 fixes angular x=-1.",
                "At lambda=0 the update is Q_j. Forced SUB at lambda=1 is identity (continuous pole extension).",
                "GREEDY freezes at lambda=1 only where Q_j is finite and its SUB branch is chosen.",
                "P / M: nine sampled metrics. Candidate is a sampled gradient; lyap/twin are finite separation growth.",
                "Quality and precision presets change grid/seed budgets. Step budget truncates the finite ensemble.",
                "Square domain: drag to pan, wheel to zoom, right-click for lambda readout. A stays in this plane.",
                "1: raw greedy. 2: raw subtract. 3: selfdual greedy. 4: compact greedy. RH host controls are retained.",
                "R: reset bounds. S: BMP. Ctrl+S/O: views. Tab: catalogue. H: help. F11: fullscreen.",
                "Finite Taylor hosts and sampled wrapper diagnostics do not certify RH or infinite-time phase boundaries."};
        }
        if(connectedQqrhLambdaMode()) {
            lines={std::string(qqrhTitle(state.generation))+" / original BKQR outer lambda wrapper",
                state.generation<=2 ? "Direct host Q=c+eta*(F(b*z+h)-F(h)); eta=0 gives the constant host c."
                                    : "Feedback host Q=z^2+c-eta*F(b*z+h)/(b*F'(b*z+h)); eta=0 gives z^2+c.",
                "Outer update: z_next=z+(1+sigma*H(lambda))*(Q-z). Inner complex eta is independent of pixel lambda.",
                "Channels: curvature, dyadic, Laplacian 1, Laplacian 2, moving boundary. Tab opens the catalogue.",
                "Curvature does not use x. Moving boundary alone uses continuous slope m; a_m=(1-m^2)/(1+m^2).",
                "B: GREEDY / SUB / ADD. C: RAW / SELFDUAL / COMPACT. K: Taylor N = 8 / 16 / 32.",
                "At lambda=0 the update is Q. Forced SUB at lambda=1 is identity by continuous pole extension.",
                "GREEDY freezes at lambda=1 only where Q is finite and its SUB branch is chosen.",
                "P / M: nine finite metrics. Candidate is a sampled gradient; lyap/twin are finite separation growth.",
                "Quality / precision choose grid and seed budgets. Steps truncate the sampled orbits.",
                "Square: drag / wheel / right-click readout. Presets 1-4 retain host controls. R: bounds. S: BMP.",
                "Ctrl+S/O: views. H: help. Formula and complete parameters are recorded with exported images.",
                "Real-variable positivity certificates do not certify complex orbits, sampled phase boundaries or RH."};
        }
        if(directQqrhMode()) {
            lines={std::string(qqrhTitle(state.generation))+" / direct iteration",
                state.generation==6 ? "K_k(y)=k*C_H(k*y)-C_H(y); z_next=c+eta*(K_k(b*z+h)-K_k(h))."
                    : state.generation<=2 ? "z_next=c+eta*(F(b*z+h)-F(h)); eta=0 gives the constant map c."
                                    : "z_next=z^2+c-eta*F(b*z+h)/(b*F'(b*z+h)); eta=0 gives z^2+c.",
                "C / B: parameter c portrait, using the fixed seed z0=0. J: Julia at selected c.",
                "Eta is an independent fixed complex coefficient. The pixel is c in the portrait and z0 in Julia.",
                "Right-click the c portrait to select c and open Julia; right-click Julia to trace a finite orbit.",
                state.generation==6 ? "C_H=-d_y[y*H'/H], H=T(-y). Kappa=1 gives exactly c; kappa=2 is dyadic. Angular x is unused."
                    : "Curvature channels omit x. Moving boundary alone uses slope m; scale-ratio curvature uses kappa.",
                "K: Taylor N = 8 / 16 / 32. 1: selected Julia. 2: c portrait. 3: N=8 Julia. 4: N=32 Julia.",
                "P: three palettes. [ / ]: hue. Comma / period: bands. Drag: pan. Wheel: zoom. R: reset.",
                "The parameter portrait samples one seed; its colours do not certify connectedness.",
                "Magnitude |z|>16 is a display cutoff. Numerical failures and finite budget limits stay unresolved.",
                "Finite-tail near-periods are recurrence observations, not certified stable cycles.",
                "S: BMP. Ctrl+S/O: views. Save orbit: finite CSV trace and parameters. HQ: 2 x 2 samples.",
                "Real-variable positivity certificates do not certify complex dynamics or RH. Tab: catalogue. H: help."};
        }
        int y = b.y + sc(9);
        for (size_t i = 0; i < lines.size() && y + sc(19) < b.y + b.h; i++, y += sc(23))
            text(dc, {b.x + sc(12), y, b.w - sc(24), sc(21)}, lines[i],
                 i == 0 ? RGB(239, 224, 255) : RGB(201, 209, 229));
    }
    void releasePaintBuffer() {
        if(paintDC && paintPrevious)SelectObject(paintDC,paintPrevious);
        if(paintBitmap)DeleteObject(paintBitmap);
        if(paintDC)DeleteDC(paintDC);
        paintDC=nullptr; paintBitmap=nullptr; paintPrevious=nullptr;
        paintWidth=paintHeight=0;
    }
    bool ensurePaintBuffer(HDC target) {
        if(paintDC && paintWidth==width && paintHeight==height)return true;
        releasePaintBuffer();
        if(width<1 || height<1)return false;
        paintDC=CreateCompatibleDC(target);
        paintBitmap=CreateCompatibleBitmap(target,width,height);
        if(!paintDC || !paintBitmap) {releasePaintBuffer();return false;}
        paintPrevious=SelectObject(paintDC,paintBitmap);
        if(!paintPrevious || paintPrevious==HGDI_ERROR) {
            paintPrevious=nullptr;releasePaintBuffer();return false;
        }
        paintWidth=width;paintHeight=height;
        return true;
    }
    void presentPaint(HDC target,HDC dc) {
        if(dc!=target)BitBlt(target,0,0,width,height,dc,0,0,SRCCOPY);
        // A UI paint timestamp, rather than GPU completion alone, includes
        // event coalescing, worker queueing and transfer back to the window.
        if(!inCatalogue && frame && frame->inputMs<0 && frame->inputAt!=Clock::time_point{}) {
            frame->inputMs=std::chrono::duration<double,std::milli>(Clock::now()-frame->inputAt).count();
            RECT status{0,canvasH,width,canvasH+sc(28)};
            InvalidateRect(window,&status,FALSE);
        }
    }
    void paint(HDC target) {
        HDC dc=ensurePaintBuffer(target)?paintDC:target;
        auto oldFont = SelectObject(dc, font);
        fill(dc, {0, 0, width, height}, RGB(6, 7, 13));
        if (inCatalogue) {
            paintCatalogue(dc);
            SelectObject(dc, oldFont);
            presentPaint(target,dc);
            return;
        }
        if (frame && !frame->pixels.empty()) {
            BITMAPINFO bi{};
            bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
            bi.bmiHeader.biWidth = frame->w;
            bi.bmiHeader.biHeight = -frame->h;
            bi.bmiHeader.biPlanes = 1;
            bi.bmiHeader.biBitCount = 32;
            bi.bmiHeader.biCompression = BI_RGB;
            SetStretchBltMode(dc, COLORONCOLOR);
            StretchDIBits(dc, 0, 0, width, canvasH, 0, 0, frame->w, frame->h, frame->pixels.data(), &bi,
                          DIB_RGB_COLORS, SRCCOPY);
        } else
            text(dc, {0, canvasH / 2 - sc(15), width, sc(30)}, "FRAKTALFAKTUR / preparing the first frame",
                 RGB(211, 192, 247), DT_CENTER);
        if ((state.family == 7 || state.family == 10) && state.plane == 0) {
            const auto &v = state.cameras[0];
            int carrierX = int(std::lround(width * (.5 - v.x / v.span)));
            if (carrierX >= 0 && carrierX < width) {
                line(dc, carrierX, 0, carrierX, canvasH, RGB(214, 191, 92));
                text(dc, {std::min(carrierX + sc(7), width - sc(300)), sc(8), sc(292), sc(19)},
                     "GRH carrier  Re(s)=1/2  <->  Re(w)=0", RGB(232, 212, 126));
            }
        }
        if (state.family == 9) {
            int pw = width / 3;
            line(dc, pw, 0, pw, canvasH, RGB(215, 205, 231), 2);
            line(dc, 2 * pw, 0, 2 * pw, canvasH, RGB(215, 205, 231), 2);
            fill(dc, {0, 0, width, sc(29)}, RGB(8, 9, 15));
            text(dc, {sc(8), sc(5), pw - sc(16), sc(20)}, "Theta_chi(z+)  phase", RGB(237, 225, 255), DT_CENTER);
            text(dc, {pw + sc(8), sc(5), pw - sc(16), sc(20)}, "R_chi(z+) = Theta'/Theta", RGB(237, 225, 255), DT_CENTER);
            text(dc, {2 * pw + sc(8), sc(5), width - 2 * pw - sc(16), sc(20)}, "literal C_chi(a,r)", RGB(237, 225, 255), DT_CENTER);
        }
        if (state.family == 11 && state.plane == 0 && state.method == 1 &&
            state.spectralD > 0 && (spectralGen2() || std::norm(state.alpha) > 0)) {
            if (spectralGen2()) {
                for (double a : {-1.0, .35, 1.0})
                    marker(dc, {0,0,width,canvasH}, state.cameras[0], C(a*state.spectralD,0), RGB(243,208,119));
            } else {
                marker(dc, {0,0,width,canvasH}, state.cameras[0], C(0,state.spectralD), RGB(243,208,119));
                marker(dc, {0,0,width,canvasH}, state.cameras[0], C(0,-state.spectralD), RGB(243,208,119));
            }
        }
        if (quantumMode()) {
            fill(dc,{0,0,width,sc(26)},RGB(8,9,17));
            text(dc,{sc(12),sc(3),width-sc(24),sc(20)},
                state.plane==0 ? "CONDITIONED HISTORIES / Bloch sphere" :
                state.plane==1 ? "EQUAL-AREA ATLAS / longitude versus Bloch height" :
                "HISTORIES + ENSEMBLE / structure retained by individual measurements",RGB(222,227,246));
        }
        if (sphereMode()) {
            fill(dc,{0,0,width,sc(26)},RGB(8,9,17));
            if (state.plane==2) {
                text(dc,{sc(12),sc(3),width/2-sc(24),sc(20)},"FRONT HEMISPHERE",RGB(222,227,246));
                text(dc,{width/2+sc(12),sc(3),width/2-sc(24),sc(20)},"BACK HEMISPHERE",RGB(222,227,246));
            } else text(dc,{sc(12),sc(3),width-sc(24),sc(20)},
                state.plane==0 ? "SPHERICAL FRACTALS / sampled limit geometry" :
                                "EQUAL-AREA ATLAS / longitude versus sphere height",RGB(222,227,246));
        }
        if (spinMode() && !spinAtlas()) {
            fill(dc,{0,0,width,sc(26)},RGB(8,9,17));
            const bool dual=state.generation!=3 && state.plane==2;
            if (dual) {
                text(dc,{sc(12),sc(3),width/2-sc(24),sc(20)},"SPIN DYNAMICS / FRONT",RGB(222,227,246));
                text(dc,{width/2+sc(12),sc(3),width/2-sc(24),sc(20)},"BACK HEMISPHERE",RGB(222,227,246));
            } else text(dc,{sc(12),sc(3),width-sc(24),sc(20)},
                (state.generation==3 ? state.plane==1 : state.plane==0)
                    ? "DISSIPATIVE SPIN / retained deterministic trajectories"
                    : "SPIN EQUAL-AREA ATLAS / longitude versus height",RGB(222,227,246));
        }
        if(cubicMode()) {
            if(state.plane==0) marker(dc,{0,0,width,canvasH},state.cameras[0],state.c,RGB(255,221,138));
            else if(diagnostic.cubicModel && diagnostic.cubicModel->validModel) {
                marker(dc,{0,0,width,canvasH},state.cameras[1],C(diagnostic.cubicModel->critical,0),RGB(255,221,138));
                marker(dc,{0,0,width,canvasH},state.cameras[1],C(-diagnostic.cubicModel->critical,0),RGB(171,151,247));
            }
        }
        if(quotientMode()) {
            if(state.plane==0) marker(dc,{0,0,width,canvasH},state.cameras[0],state.c,RGB(255,221,138));
            else if(diagnostic.quotientModel && diagnostic.quotientModel->validModel) {
                marker(dc,{0,0,width,canvasH},state.cameras[1],C(0,1),RGB(255,221,138));
                marker(dc,{0,0,width,canvasH},state.cameras[1],C(0,-1),RGB(171,151,247));
                marker(dc,{0,0,width,canvasH},state.cameras[1],C(0,0),RGB(255,120,120));
            }
        }
        if(tripleMode()) {
            if(state.plane==0)marker(dc,{0,0,width,canvasH},state.cameras[0],state.c,RGB(255,221,138));
            else if(diagnostic.tripleModel&&diagnostic.tripleModel->validModel) {
                const auto starts=diagnostic.tripleModel->criticalPoints(state.c);
                const COLORREF colours[]={RGB(255,221,138),RGB(171,151,247),RGB(114,218,218)};
                for(size_t k=0;k<starts.size();++k)marker(dc,{0,0,width,canvasH},state.cameras[1],starts[k],colours[k]);
                if(state.generation==3)marker(dc,{0,0,width,canvasH},state.cameras[1],C{},RGB(255,120,120));
            }
        }
        if(directQqrhMode() && state.plane==0)
            marker(dc,{0,0,width,canvasH},state.cameras[0],state.c,RGB(255,221,138));
        if(masterMode()) {
            if(state.plane==0)
                marker(dc,{0,0,width,canvasH},state.cameras[0],state.c,RGB(255,221,138));
            else if(state.plane==2) {
                const auto &cam=state.cameras[2];
                const double axisX=width*(.5-cam.x/cam.span);
                const double axisY=canvasH*.5+cam.y*width/cam.span;
                if(std::isfinite(axisX) && axisX>=0 && axisX<width)
                    line(dc,int(axisX),0,int(axisX),canvasH,RGB(86,91,116));
                if(std::isfinite(axisY) && axisY>=0 && axisY<canvasH)
                    line(dc,0,int(axisY),width,int(axisY),RGB(86,91,116));
                text(dc,{sc(12),sc(6),sc(260),sc(20)},"Im lambda / fixed seed z0=0",RGB(216,218,235));
                text(dc,{width-sc(130),canvasH-sc(27),sc(118),sc(20)},"Re lambda",RGB(216,218,235),DT_RIGHT);
                marker(dc,{0,0,width,canvasH},cam,state.alpha,RGB(255,221,138));
                marker(dc,{0,0,width,canvasH},cam,C(1,0),RGB(174,163,252));
                const double refX=width*(.5+(1-cam.x)/cam.span);
                if(std::isfinite(refX) && std::isfinite(axisY) && refX>=7 && refX<width-7 && axisY>=7 && axisY<canvasH-7)
                    text(dc,{std::clamp(int(refX)+sc(9),sc(8),std::max(sc(8),width-sc(245))),
                             std::clamp(int(axisY)-sc(26),sc(5),std::max(sc(5),canvasH-sc(25))),sc(235),sc(20)},
                         "lambda=1 / reference only",RGB(197,187,255));
            }
        }
        drawLegacyLambdaGuides(dc);
        drawOrbit(dc);
        hilbertInspector(dc);
        helpOverlay(dc);
        fill(dc, {0, canvasH, width, height - canvasH}, RGB(16, 17, 25));
        line(dc, 0, canvasH, width, canvasH, RGB(98, 78, 125));
        std::string plane = state.family == 9 ? "shared (log2 a, logit r)" : state.plane == 0   ? "z-start"
                            : state.plane == 1 ? (state.family == 1   ? "a plane"
                                                  : state.family == 4 ? "feedback plane"
                                                  : state.family == 3 ? "q embedding plane"
                                                  : state.family >= 2 ? "target plane"
                                                                      : "c plane")
                                               : (state.family == 1   ? "b plane"
                                                  : state.family >= 2 ? "relaxation plane"
                                                                      : "alpha plane");
        if(cubicMode()) plane=state.plane==0 ? "critical-orbit parameter atlas" : "selected cubic Julia plane";
        if(tripleMode()) plane=state.plane==0 ? "all-critical parameter atlas" : state.generation==2 ? "return-map Julia plane" : "selected Julia / infinity basin";
        if(quotientMode()) plane=state.plane==0 ? "critical escape parameter atlas" : "rational Julia / infinity basin";
        if(masterMode()) plane=state.plane==0 ? "c portrait / fixed seed z0=0" : state.plane==2 ?
            "complex lambda / fixed seed z0=0" : "selected Julia / finite orbit diagnostics";
        if(legacyLambdaMode())plane="original BKQR lambda plane / square domain";
        if(connectedLambdaMode())plane="BKQR outer lambda / RH host / square domain";
        if(connectedQqrhLambdaMode())plane="BKQR outer lambda / QQRH host / square domain";
        if(directQqrhMode())plane=state.plane==0 ? "c portrait / fixed seed z0=0 / fixed eta" : "selected Julia / direct QQRH";
        if(extraMode()) {
            const int driver=state.generation+17;
            plane=frk_extra::is_field_driver(driver) ? "complex z plane / direct Hilbert geometry" :
                frk_extra::seed_at_pixel(driver) ? "complex starting-z plane" :
                state.plane==0 ? "parameter c plane / seed z0=0" : "Julia z plane / selected c";
        }
        if(projectMode() || riemannMode()) plane="complex initial-state plane";
        if (state.family == 11)
            plane = state.plane == 0 ? "x-start plane" : state.plane == 1 ? "c / selected-seed slice" : "eta / selected-seed slice";
        if (quantumMode()) plane=state.plane==0 ? "Bloch sphere" : state.plane==1 ? "Equal-area atlas" : "Histories / ensemble";
        if (spinMode()) plane=spinAtlas() ? "Twist / damping atlas" : state.generation==3 ? (state.plane==1 ? "Selected spin sphere" : "Selected equal-area atlas") : state.plane==0 ? "Spin sphere" : state.plane==1 ? "Equal-area atlas" : "Front / back hemispheres";
        if (sphereMode()) plane=state.plane==0 ? "Riemann sphere" : state.plane==1 ? "Equal-area atlas" : "Front / back hemispheres";
        std::string quality =
            frame && frame->generation == generation.load() &&
                (!masterMode() || samePaintScene(frame->state,state))
                ? (frame->stride == 1 ? "FULL RES" : std::to_string(frame->stride) + "x PREVIEW")
                : "REFINING";
        std::string top = familyName(state.family) +
                          (state.family ? " g" + std::to_string(state.generation) : "") + "  |  " + plane +
                          "  |  " + quality + (masterMode() && frame && frame->gpuMetrics.cacheHit ? " / cached colour" : "") + "  |  " + (frame ? num(frame->ms, 1) : "--") + " ms" +
                          (masterMode() && frame && frame->inputMs>=0 ? " / input "+num(frame->inputMs,0)+" ms" : "") + "  |  " +
                          backend;
        text(dc, {sc(12), cy(4), width - sc(185), sc(22)}, top, RGB(235, 229, 248));
        button(dc, {width - sc(165), cy(4), sc(153), sc(23)}, "Renderers [Tab]");
        std::string info =
            (legacyLambdaMode() || directAtlasMode() || projectMode() || riemannMode() || state.family == 11 || sampleMode()) ? std::string{} : state.family ? "Generation " + std::to_string(state.generation) + "    degree " +
                               std::to_string(diagnostic.gardenData->degree) + "    " +
                               (state.method == 1 ? "Halley" : "Newton") + "    " +
                               (state.family == 1   ? "a = "
                                : state.family == 4 ? "feedback = "
                                : state.family == 3 ? "q = "
                                                    : "target = ") +
                               complexText(state.c) + "    " +
                               (state.family == 1 ? "b = " : "relaxation = ") + complexText(state.alpha)
                         : "c = " + complexText(state.c) + "    alpha = " + complexText(state.alpha) +
                               "    critical periods: + " + periodText(diagnostic.a.period) + " / - " +
                               periodText(diagnostic.b.period) + (diagnostic.same ? " (same cycle)" : "");
        if (state.family == 6)
            info = lfunctionName(state.generation) + "    " + (state.method ? "Halley" : "Newton") +
                   "    target = " + complexText(state.c) + "    relaxation = " + complexText(state.alpha);
        if (state.family == 7)
            info = "primitive chi mod 5    Cayley sheet k=" + std::to_string(state.generation) + "    " +
                   (state.method ? "Halley" : "Newton") + "    normalized target = " + complexText(state.c) +
                   "    relaxation = " + complexText(state.alpha);
        if (state.family == 8)
            info = "primitive chi5 x chibar5 pair    " + std::string(state.method ? "Halley" : "Newton") +
                   "    normalized target = " + complexText(state.c) + "    relaxation = " + complexText(state.alpha);
        if (state.family == 9)
            info = "primitive chi5 completed-pair closure    x=log2(a)    y=logit(r)    three synchronized exact-observable panels";
        if (state.family == 10)
            info = "primitive chi mod 5    Cayley sheet 3    " + std::string(state.method ? "Halley" : "Newton") +
                   " + phase feedback    kappa=" + num(state.c.real(),3) + "    tau=" + num(state.c.imag(),3);
        if (state.family == 11)
            info = "Direct map   c=" + complexText(state.c) + "   eta=" + complexText(state.alpha) +
                   "   d=" + num(state.spectralD,3) + "   seed=" + complexText(state.seed);
        if (spectralGen2())
            info = "Three atoms / jet " + std::to_string(state.spectralOrder) + "   c=" + complexText(state.c) +
                   "   eta=" + complexText(state.alpha) + "   d=" + num(state.spectralD,3) +
                   "   middle weight=" + num(state.spectralWeight,3);
        if (quantumMode()) {
            info=std::string(quantum::geometryName(state.quantum.geometry))+" / strength="+num(state.quantum.strength,3)+
                " / kick="+num(state.quantum.kick,3)+" / replay seed="+std::to_string(state.quantum.randomSeed);
            if (frame && frame->generation==generation.load()) info+=" / samples="+std::to_string(frame->quantumStats.samples)+
                " / |history mean|="+num(frame->quantumStats.meanNorm,4);
        }
        if (sphereMode()) {
            info=state.generation==3 ? "McMullen / log10|lambda|="+num(state.sphere.lambdaLog,3)+" / arg="+num(state.sphere.lambdaPhase,3)
                 : state.generation==2 ? "Schottky / circle gap="+num(state.sphere.gap,4) : "Apollonian / tangent reflection circles";
            info+=" / display boost="+num(state.sphere.boost,3)+" / seed="+std::to_string(state.quantum.randomSeed);
            if (frame && frame->generation==generation.load()) info+=" / samples="+std::to_string(frame->sphereStats.samples);
        }
        if (spinMode()) {
            info="k="+num(state.spin.twist,3)+" / beta="+num(state.spin.tilt,3)+" / tau="+num(state.spin.damping,3);
            if (frame && frame->generation==generation.load()) {
                const auto &st=frame->spinStats;
                info+=spinAtlas() ? " / grid="+std::to_string(st.gridWidth)+"x"+std::to_string(st.gridHeight)+
                                   " / starts/cell="+std::to_string(st.startsPerCell)
                                 : " / samples="+std::to_string(st.samples);
                info+=" / finite-time growth="+num(st.meanLyapunov,4)+" / log area="+num(st.meanAreaRate,4);
            }
        }
        if(projectMode()) info=std::string(projectTitle(state.generation))+" / "+
            (state.generation==1 ? "finite approximation degree "+std::to_string(state.project.grhDegree)
             : state.generation==2 ? std::to_string(state.project.piTerms)+" actual pi pairs / offset="+std::to_string(state.project.piOffset)
             : state.generation==3 ? "nu*t="+num(state.project.nseTime,6)+" / modes="+std::to_string(state.project.nseModes)
             : state.generation==4 ? "F composed with F / finite orbit diagnostics"
             : "q=2 / p="+num(state.project.pcdProbability,9));
        if(projectMode() && diagnostic.projectModel && !diagnostic.projectModel->validModel)
            info="MODEL UNAVAILABLE: "+diagnostic.projectModel->statusText;
        if(riemannMode()) {
            info=(state.generation>=4 ? std::string("Five fixed H2 moments / ")
                                     : "Finite H"+std::to_string(state.riemann.order)+" witness / ")+
                std::string(state.method ? "Halley" : "Newton")+" / relaxation="+complexText(state.alpha);
            if(diagnostic.riemannModel) {
                info+=" / degree="+std::to_string(diagnostic.riemannModel->degree);
                if(!diagnostic.riemannModel->validModel) info="MODEL UNAVAILABLE: "+diagnostic.riemannModel->statusText;
            }
        }
        if(cubicMode()) {
            info="Direct cubic / c="+complexText(state.c);
            if(diagnostic.cubicModel) {
                info+=" / A="+num(diagnostic.cubicModel->A,8)+" / critical points +/- "+num(diagnostic.cubicModel->critical,6);
                if(!diagnostic.cubicModel->validModel) info="MODEL UNAVAILABLE: "+diagnostic.cubicModel->statusText;
            } else info="MODEL UNAVAILABLE: missing cubic model";
        }
        if(quotientMode()) {
            info="Jacobi quotient / c="+complexText(state.c);
            if(diagnostic.quotientModel) {
                info+=" / lambda="+num(diagnostic.quotientModel->lambda,7)+" / critical +/-i; pole 0";
                if(!diagnostic.quotientModel->validModel) info="MODEL UNAVAILABLE: "+diagnostic.quotientModel->statusText;
            } else info="MODEL UNAVAILABLE: missing quotient model";
        }
        if(extraMode()) {
            const auto& entry=frk_extra::extra_catalogue[state.generation-1];
            info=std::string(entry.label)+" / real lambda="+num(state.alpha.real(),4)+
                " / selected c="+complexText(state.c)+" / "+
                (state.plane==0?"parameter c":"starting z")+
                " / finite orbit classification";
            if(hilbertMode()) info=std::string(entry.label)+" / N=12 / R="+num(state.alpha.real(),4)+
                " / physical r="+num(200.*state.alpha.real(),2)+
                (state.generation==9?" / squared Hilbert overlap [0,1]":
                 state.generation==12?" / projective metric color g/(1+g)":" / finite Hilbert dynamics");
        }
        if(tripleMode()) {
            info=std::string(momenttriple::title(state.generation))+" / c="+complexText(state.c);
            if(diagnostic.tripleModel) {
                if(!diagnostic.tripleModel->validModel)info="MODEL UNAVAILABLE: "+diagnostic.tripleModel->statusText;
                else info+=state.generation==2 ? " / one step = two phases" : " / all finite critical orbits";
            } else info="MODEL UNAVAILABLE: missing three-way model";
        }
        if(masterMode()) {
            info=std::string(masterTitle(state.generation))+" / c="+complexText(state.c)+
                " / N="+std::to_string(8<<(state.spectralOrder-1))+" / lambda="+complexText(state.alpha)+
                " / "+rhmaster::mapName(state.rhExperiment.map)+" / "+rhmaster::chartName(state.rhExperiment.chart)+
                " / "+rhmaster::colouringName(state.rhExperiment.coloring)+
                " / R="+num(state.seed.real(),3)+" / level="+num(value(9),3)+
                " / h="+num(state.spectralD,3)+" / b="+num(state.seed.imag(),3);
            if(!diagnostic.masterModel || !diagnostic.masterModel->validModel)
                info="MODEL UNAVAILABLE: "+(diagnostic.masterModel ? diagnostic.masterModel->statusText : std::string("missing candidate model"));
        }
        if(directQqrhMode()) {
            info=std::string(qqrhTitle(state.generation))+" / direct / c="+complexText(state.c)+
                " / eta="+complexText(state.alpha)+" / N="+std::to_string(8<<state.method)+
                " / x="+(state.generation==1 || state.generation==6 ? std::string("unused") : num(state.seed.real(),3))+
                " / h="+num(state.spectralD,3)+" / b="+num(state.seed.imag(),3);
            if(state.generation==5)info+=" / m="+num(state.qqrhSlope,3);
            if(state.generation==6)info+=" / kappa="+num(state.qqrhScale,3);
            if(!diagnostic.qqrhDirectModel || !diagnostic.qqrhDirectModel->validModel)
                info="MODEL UNAVAILABLE: "+(diagnostic.qqrhDirectModel ? diagnostic.qqrhDirectModel->statusText : std::string("missing direct QQRH model"));
        }
        if(legacyLambdaMode()) {
            const auto settings=bkqrSettings(state);
            info=std::string(bkqrlambda::hostName(settings.host))+" / "+bkqrlambda::branchName(settings.branch)+
                " / "+bkqrlambda::couplingName(settings.coupling)+" / c="+complexText(state.c)+
                " / power="+std::to_string(settings.power)+" / metric="+legacyMetricTitle(settings.metric)+
                " / grid="+std::to_string(bkqrlambda::gridSize(settings))+" / seeds="+std::to_string(bkqrlambda::seedCount(settings));
        }
        if(connectedLambdaMode()) {
            const auto settings=bkqrSettings(state);
            info=std::string(connectedQqrhLambdaMode() ? qqrhTitle(state.generation) : masterTitle(state.generation))+" / inner eta="+complexText(state.alpha)+
                " / N="+std::to_string(8<<state.method)+" / "+bkqrlambda::branchName(state.orbitMode)+
                " / H="+bkqrlambda::couplingName(state.spectralOrder-1)+" / metric="+legacyMetricTitle(state.palette)+
                " / grid="+std::to_string(bkqrlambda::gridSize(settings))+" / seeds="+std::to_string(bkqrlambda::seedCount(settings));
            if(connectedQqrhLambdaMode() && state.generation==5)info+=" / m="+num(state.qqrhSlope,3);
        }
        text(dc, {sc(12), cy(28), width - sc(24), sc(17)}, info, RGB(165, 184, 209));
        static const char *ampPresets[] = {"1  Fourfold", "2  Seventh", "3  Triplet", "4  Frontier"};
        static const char *gpPresets[] = {"1  Balanced", "2  Alternating", "3  Twisted", "4  Quartic branch"};
        static const char *percoPresets[] = {"1  Target 1/2", "2  Target 1/4", "3  Target 3/4",
                                             "4  Complex target"};
        static const char *crownPresets[] = {"1  Binary crown", "2  Spiral digits", "3  Folded crown",
                                             "4  Twisted step"};
        static const char *restartPresets[] = {"1  Complex feedback", "2  Markov return", "3  Phase reversal",
                                               "4  Twisted step"};
        static const char *cayleyPresets[] = {"1  Complex target", "2  Half transmission", "3  Halley branches", "4  Twisted step"};
        static const char *dirichletPresets[] = {"1  Newton zeros", "2  Halley zeros", "3  Complex level", "4  Twisted step"};
        static const char *grhCayleyPresets[] = {"1  Primitive zeros", "2  Halley basins", "3  Boundary zoom", "4  Twisted step"};
        static const char *pairPresets[] = {"1  Paired zeros", "2  Halley mirror", "3  Central weave", "4  Twisted pair"};
        static const char *observablePresets[] = {"1  Boundary zeros", "2  Wide domain", "3  Higher zeros", "4  Extreme boundary"};
        static const char *feedbackPresets[] = {"1  Flowing basins", "2  Eddy zoom", "3  Turbulent weave", "4  Halley feedback"};
        static const char *spectralPresets[] = {"1  Spectral lace", "2  Islands", "3  Filaments", "4  Quadratic control"};
        static const char *spectralGen2Presets[] = {"1  Triad", "2  Split branches", "3  Second jet", "4  Third jet"};
        static const char *quantumPresets[]={"1  Tetrahedron", "2  Octahedron", "3  Icosahedron", "4  Kicked dodecahedron"};
        static const char *apollonianPresets[]={"1  Symmetric gasket", "2  Conformal squeeze", "3  Dual hemispheres", "4  Gasket atlas"};
        static const char *schottkyPresets[]={"1  Dense dust", "2  Separated islands", "3  Sparse dust", "4  Tangency limit"};
        static const char *mcmullenPresets[]={"1  Published rings", "2  Fine rings", "3  Complex lambda", "4  Larger lambda"};
        static const char *spinPresets[]={"1  Published attractor", "2  Weak damping", "3  Strong twist", "4  Conservative control"};
        static const char *stretchPresets[]={"1  Published stretching", "2  Strong twist", "3  Damping control", "4  Area-preserving control"};
        static const char *spinAtlasPresets[]={"1  Broad atlas", "2  Fine structure", "3  Tilted atlas", "4  Strong twist"};
        static const char *triplePresets[]={"1  Critical atlas","2  c=0 Julia","3  Critical-fixed Julia","4  Complex Julia"};
        static const char *masterPresets[]={"1  Default Julia","2  Seeded parameter","3  Degree 8 Julia","4  Degree 32 Julia"};
        std::array<std::string,4> labPresetLabels;
        std::array<const char *,4> labPresets;
        for(int i=0;i<4;++i) {
            labPresetLabels[i]=std::to_string(i+1)+"  "+focus_presets::presetName(focusBank,i);
            labPresets[i]=labPresetLabels[i].c_str();
        }
        static const char *directQqrhPresets[]={"1  Selected Julia","2  Seeded parameter","3  Degree 8 Julia","4  Degree 32 Julia"};
        static const char *legacyPresets[]={"1  Raw greedy","2  Raw subtract","3  Selfdual greedy","4  Compact greedy"};
        const bool fixedExtra=extraMode() && frk_extra::seed_at_pixel(state.generation+17);
        const char *extraPresets[]={"1  Standard z-plane", fixedExtra ? "2  Reset parameters" : "2  Parameter c-plane",
            fixedExtra ? "3  Current parameters" : "3  Selected Julia/c", "4  Shifted palette"};
        if(hilbertMode()) {
            extraPresets[0]="1  Default R=0.65";
            if(hilbertField()) {
                extraPresets[1]="2  Radius R=0.25";extraPresets[2]="3  Radius R=1.10";
                extraPresets[3]="4  z=1 close-up";
            } else if(fixedExtra) {
                extraPresets[1]="2  Radius R=0.35";extraPresets[2]="3  Radius R=1.10";
                extraPresets[3]="4  Boundary close-up";
            } else extraPresets[3]="4  Complex Julia R=0.85";
        }
        static const char *cubicPresets[]={"1  Critical atlas", "2  Symmetric Julia", "3  Critical-fixed Julia", "4  Complex Julia"};
        static const char *riemannPresets[7][4]={
            {"1  Three-node witness", "2  Node close-up", "3  Two-node H1", "4  Complex relaxation"},
            {"1  All cubic branches", "2  Canonical small branch", "3  Fold near x=1", "4  Halley branches"},
            {"1  Reciprocal atlas", "2  Inner cubic branches", "3  Outer sheets", "4  Complex relaxation"},
            {"1  Complex pencil", "2  Central detail", "3  Halley pencil", "4  Complex relaxation"},
            {"1  Mirrored kernel", "2  Central detail", "3  Halley kernel", "4  Complex relaxation"},
            {"1  Jacobi closure", "2  Central detail", "3  Halley spectrum", "4  Complex relaxation"},
            {"1  Joukowski sheets", "2  Unit-circle detail", "3  Halley sheets", "4  Complex relaxation"}};
        static const char *projectPresets[5][4]={
            {"1  Pole field", "2  Complex feedback", "3  Wider carrier", "4  Fine carrier"},
            {"1  64 terms", "2  32 terms", "3  128 terms", "4  Shifted block"},
            {"1  Six modes", "2  Young spectrum", "3  Damped spectrum", "4  Steeper spectrum"},
            {"1  Cycle basins", "2  Wide return field", "3  Left boundary", "4  Fine return field"},
            {"1  d=3 anchor", "2  d=4 anchor", "3  d=6 anchor", "4  Complex branches"}};
        for (int i = 0; i < 4; i++)
            button(dc, row(i, 4, 49),
                   (legacyLambdaMode() ? legacyPresets : directQqrhMode() ? directQqrhPresets : masterMode() ? (labControls ? labPresets.data() : masterPresets) : tripleMode() ? triplePresets : extraMode() ? extraPresets : directAtlasMode() ? cubicPresets : riemannMode() ? riemannPresets[std::clamp(state.generation,1,7)-1] : projectMode() ? projectPresets[std::clamp(state.generation,1,5)-1] : spinMode() ? (state.generation==1 ? spinPresets : state.generation==2 ? stretchPresets : spinAtlasPresets)
                    : sphereMode() ? (state.generation==1 ? apollonianPresets : state.generation==2 ? schottkyPresets : mcmullenPresets)
                    : quantumMode() ? quantumPresets : spectralGen2() ? spectralGen2Presets : state.family == 11 ? spectralPresets
                    : state.family == 10  ? feedbackPresets
                    : state.family == 9   ? observablePresets
                    : state.family == 8 ? pairPresets
                    : state.family == 7 ? grhCayleyPresets
                    : state.family == 6 ? dirichletPresets
                    : state.family == 5 ? cayleyPresets
                    : state.family == 1 ? gpPresets
                    : state.family == 4 ? restartPresets
                    : state.family == 3 ? crownPresets
                    : state.family == 2 ? percoPresets
                                        : ampPresets)[i],
                   masterMode() && labControls ? labPresetSelected(i) : state.preset == i);
        if (!state.family) {
            static const char *modes[] = {"Plane: z-start", "Plane: c",   "Plane: alpha", "Both critical",
                                          "+ critical",     "- critical", "Seed"};
            for (int i = 0; i < 7; i++)
                button(dc, row(i, 7, 81), modes[i], i < 3 ? state.plane == i : state.orbitMode == i - 3,
                       i < 3 || state.plane != 0);
        } else if(directQqrhMode()) {
            const char *modes[]={"Seeded c plane [C/B]","Julia plane [J]","Cycle Taylor N [K]",
                qqrhTitle(state.generation),"Direct Q iteration","Formula / Help"};
            for(int i=0;i<6;++i)button(dc,row(i,6,81),modes[i],(i<2 && state.plane==i) || (i==5 && help),i<3 || i==5);
        } else if(connectedLambdaMode()) {
            const char *modes[]={"Outer lambda plane [A]",connectedQqrhLambdaMode() ? qqrhTitle(state.generation) : state.generation==1 ? "Master RH host" : state.generation==2 ? "Sub_1 RH host" : "Sub_2 RH host",
                "Branch [B]","H coordinates [C]","Cycle Taylor N [K]","Formula / Help"};
            for(int i=0;i<6;++i)button(dc,row(i,6,81),modes[i],i==0 || (i==5 && help),i!=1);
        } else if(legacyLambdaMode()) {
            const char *modes[]={"Lambda plane [A]","Cycle host [F]","Branch [B]","Coordinates [C]","Metric [P/M]","Formula / Help"};
            for(int i=0;i<6;++i)button(dc,row(i,6,81),modes[i],i==0 || (i==5 && help));
        } else if(masterMode()) {
            if(labControls) {
                const std::string mapLabel=std::string(rhmaster::mapName(state.rhExperiment.map))+" [F]";
                const std::string chartLabel=std::string(rhmaster::chartName(state.rhExperiment.chart))+" [D]";
                const std::string colorLabel=std::string(rhmaster::colouringName(state.rhExperiment.coloring))+" [I]";
                const std::string bankLabel=std::string(focus_presets::bankName(focusBank))+" [T]";
                const char *modes[]={mapLabel.c_str(),chartLabel.c_str(),colorLabel.c_str(),
                    bankLabel.c_str(),"Main controls [V]","Formula / Help"};
                for(int i=0;i<6;++i)button(dc,row(i,6,81),modes[i],i==4 || (i==5 && help));
            } else {
                const char *modes[]={"Seeded c plane [C/B]", "Julia plane [J]", "Lambda plane [A]",
                    "Cycle degree [K]", "Lab controls [V]", "Formula / Help"};
                for(int i=0;i<6;++i)button(dc,row(i,6,81),modes[i],(i<3 && state.plane==i) || (i==5 && help));
            }
        } else if(extraMode()) {
            const bool seed=frk_extra::seed_at_pixel(state.generation+17);
            const char *modes[]={"Parameter c-plane [C/B]", "Starting z-plane [J]",
                hilbertField() ? (state.generation==12 ? "Projective metric g" : "Hilbert overlap [0,1]") : seed ? "Fixed root iteration" : "Finite escape portrait",
                frk_extra::extra_catalogue[state.generation-1].label,
                hilbertMode() ? "Finite Hardy model N=12" : "Direct numerical dynamics", "Formula / Help"};
            for(int i=0;i<6;++i)
                button(dc,row(i,6,81),modes[i],(i<2&&state.plane==i)||(i==5&&help),
                       (i==0&&!seed)||i==1||i==5);
        } else if(directAtlasMode()) {
            const char *modes[]={"Critical atlas [C/B]", "Julia plane [J]", "Both critical orbits",
                tripleMode() ? "Direct moment iteration" : quotientMode() ? "Direct rational iteration" : "Direct cubic iteration", "Fixed RH moments", "Formula / Help"};
            for(int i=0;i<6;++i) button(dc,row(i,6,81),modes[i],(i<2 && state.plane==i) || (i==5 && help),i<2 || i==5);
        } else if(riemannMode()) {
            const char *modes[]={"Initial-state plane [J]", state.generation>=4 ? "Five fixed H2 moments" : "Finite witness nodes",
                state.generation>=4 ? "Fixed moment polynomial" : "Polynomial pullback",
                "Method: Newton", "Method: Halley", "Formula / Help"};
            for(int i=0;i<6;++i) button(dc,row(i,6,81),modes[i],i==0 || (i==3 && state.method==0) ||
                (i==4 && state.method==1) || (i==5 && help),i==0 || i>=3);
        } else if(projectMode()) {
            const char *modes[]={"Initial-state plane [J]", "Fixed update rule", "No parameter plane",
                state.generation==1 || state.generation==4 ? "Direct iteration" : "Newton iteration", "Project controls below", "Formula / Help"};
            for(int i=0;i<6;++i) button(dc,row(i,6,81),modes[i],i==0 || (i==5 && help),i==0 || i==5);
        } else if (spinMode()) {
            const char *modes[]={state.generation==3 ? "Parameter atlas [J]" : "Sphere [J]",
                state.generation==3 ? "Selected sphere [C]" : "Equal-area atlas [C]",
                state.generation==3 ? "Selected atlas [A]" : "Front / back [A]",
                spinAtlas() ? "Atlas guides" : "Sphere guides", "New initial points [N]", "Formula / Help"};
            for(int i=0;i<6;++i) button(dc,row(i,6,81),modes[i],i<3 ? state.plane==i : i==3 ? state.method==1 : i==5 && help);
        } else if (sphereMode()) {
            static const char *modes[]={"Sphere [J]", "Equal-area atlas [C]", "Front / back [A]", "Sphere guides", "New sample [N]", "Formula / Help"};
            for(int i=0;i<6;++i) button(dc,row(i,6,81),modes[i],i<3 ? state.plane==i : i==3 ? state.method==1 : i==5 && help);
        } else if (quantumMode()) {
            static const char *modes[]={"Sphere [J]", "Equal-area atlas [C]", "Mean / histories [A]", "Sphere guides", "New history [N]", "Formula / Help"};
            for(int i=0;i<6;++i) button(dc,row(i,6,81),modes[i],i<3 ? state.plane==i : i==3 ? state.method==1 : i==5 && help);
        } else if (state.family == 11) {
            std::array<std::string,6> modes = {"Plane: x-start", "Plane: c", "Plane: eta",
                                              spectralGen2() ? "Atom guides" : "Pole markers", "Seed = 0", "Formula / Help"};
            for (int i=0;i<6;i++)
                button(dc,row(i,6,81),modes[i], i<3 ? state.plane==i : i==3 ? state.method==1 : i==5 && help);
        } else if (state.family == 9) {
            std::array<std::string, 6> modes = {"Shared theorem plane", "LEFT: Theta phase", "MID: log derivative",
                                                "RIGHT: observable", "No Newton step", "Formula / Help"};
            for (int i = 0; i < 6; i++)
                button(dc, row(i, 6, 81), modes[i], i == 0 || (i == 5 && help), i == 0 || i == 5);
        } else {
            std::array<std::string, 6> modes = {"Plane: z-start",
                                                state.family == 1   ? "Plane: a"
                                                : state.family == 4 ? "Plane: feedback"
                                                : state.family == 3 ? "Plane: q"
                                                : state.family == 10 ? "Plane: phase feedback"
                                                                    : "Plane: target",
                                                state.family == 1 ? "Plane: b" : "Plane: relaxation",
                                                "Method: Newton",
                                                "Method: Halley",
                                                "Formula / Help"};
            for (int i = 0; i < 6; i++)
                button(dc, row(i, 6, 81), modes[i],
                       i < 3   ? state.plane == i
                       : i < 5 ? state.method == i - 3
                               : help);
        }
        static const char *ampPalettes[] = {"Basins / period bands", "Mean stretch", "Pole orbit trap"};
        static const char *gardenPalettes[] = {"Root phase / basins", "Iteration bands", "Relative residual"};
        static const char *crownPalettes[] = {"Cylinder basins", "Iteration bands", "Nearest-root distance"};
        static const char *dirichletPalettes[] = {"Root / level basins", "Iteration bands", "Absolute residual"};
        static const char *observablePalettes[] = {"Balanced contrast", "Boundary boost", "Signed-field boost"};
        static const char *feedbackPalettes[] = {"Flow basins / ribbons", "Iteration bands", "Residual structure"};
        static const char *spectralPalettes[] = {"Smooth escape", "Pole orbit trap", "Orbit phase"};
        static const char *spectralGen2Palettes[] = {"Smooth escape", "Denominator trap", "Orbit phase"};
        static const char *quantumPalettes[]={"Visit density", "Measurement outcome", "Bloch height"};
        const char *spherePalettes[]={"Visit density", state.generation==3 ? "Inverse branch" : "Generating circle", "Sphere height"};
        const char *spinPalettes[]={spinAtlas() ? "Finite-time Lyapunov" : "Visit density",
            spinAtlas() ? "Mean log area change" : "Finite-time Lyapunov", spinAtlas() ? "Final spin height" : "Mean log area change"};
        const char *extraPalettes[]={"Basin / escape phase", "Shifted orbit colors", "Reversed orbit colors"};
        const char *cubicPalettes[]={state.plane==0 ? "Critical fates" : "Escape phase", "Escape time", "Finite orbit structure"};
        const char *masterPalettes[]={state.plane==2 ? "Finite period / orbit phase" : "Orbit phase", "Iteration bands", "Finite orbit structure"};
        const char *projectPalettes[]={"Basin / cycle phase", "Convergence / recurrence time", "Finite orbit structure"};
        const char *legacyMetrics[]={"Previous metric",legacyMetricTitle(state.palette),"Next metric [P/M]"};
        for (int i = 0; i < 3; i++)
            button(dc, row(i, 3, 113),
                   (legacyLambdaMode() ? legacyMetrics : (masterMode() || directQqrhMode()) ? masterPalettes : extraMode() ? extraPalettes : directAtlasMode() ? cubicPalettes : (projectMode() || riemannMode()) ? projectPalettes : spinMode() ? spinPalettes : sphereMode() ? spherePalettes : quantumMode() ? quantumPalettes : spectralGen2() ? spectralGen2Palettes : state.family == 11 ? spectralPalettes
                    : state.family == 10 ? feedbackPalettes
                    : state.family == 9 ? observablePalettes
                    : (state.family == 6 || state.family == 7 || state.family == 8) ? dirichletPalettes
                    : state.family == 3 ? crownPalettes
                    : state.family    ? gardenPalettes
                                      : ampPalettes)[i],
                   legacyLambdaMode() ? i==1 : state.palette == i);
        const char *actions[] = {spinAtlas() ? "Reset atlas bounds" : "Reset camera",
                                 "Reset all",
                                 "Save BMP",
                                 "Open renders",
                                 sampleMode() ? "Save stats" : "Save orbit",
                                 spinMode() ? "New initial points" : sphereMode() ? "New sample" : quantumMode() ? "New history" : "Clear orbit",
                                 help ? "Hide help" : "Help",
                                 "Quit"};
        for (int i = 0; i < 8; i++)
            button(dc, row(i, 8, 146, 29), actions[i], false,
                   i == 4 || i == 5 ? (sampleMode() ? (i==5 || (frame && frame->stride==1 && frame->generation==generation.load())) : orbitVisible)
                   : i == 2         ? frame && frame->stride == 1 && frame->generation == generation.load()
                                    : true);
        const bool spinScan=spinMode() && state.generation==3;
        const char *explore[] = {legacyLambdaMode() ? "Lambda plane only" : directAtlasMode() ? "Selected Julia [E/J]" : spinScan ? "Explore selected [E]" : sampleMode() ? "More samples" : "Explore point",
                                 legacyLambdaMode() ? "Fixed finite ensemble" : (masterMode() || directQqrhMode()) ? "Seeded parameter [B]" : (directAtlasMode() || spinScan) ? "Parameter c-plane [B]" : sampleMode() ? "Fewer samples" : "Back", "Save view", "Load view",
                                 hqBusy ? "HQ rendering..." : legacyLambdaMode() ? "Save large BMP" : spinAtlas() ? "HQ (denser grid)" : sampleMode() ? "HQ (more samples)" : "Save HQ (2 x 2)"};
        for (int i = 0; i < 5; i++)
            button(dc, row(i, 5, 182, 29), explore[i], false,
                   i == 0   ? (!legacyLambdaMode() && ((directAtlasMode() || spinScan) ? state.plane!=1 : sampleMode() ? state.quantum.samplePower<23 : state.plane != 0))
                   : i == 1 ? (!legacyLambdaMode() && (!(extraMode() && frk_extra::seed_at_pixel(state.generation+17))) && ((directAtlasMode() || spinScan) ? state.plane!=0 : sampleMode() ? state.quantum.samplePower>14 : backValid))
                   : i == 4 ? !hqBusy
                            : true);
        for (auto s : sliders())
            slider(dc, s);
        if(connectedRhLambdaMode()) {
            const auto b=row(2,3,358,38);
            text(dc,b,"Outer lambda = pixel; eta = host",RGB(165,184,209),DT_CENTER);
        }
        if(directQqrhMode()) {
            const auto b=row(1,3,358,38);
            text(dc,{b.x,b.y,width-b.x-sc(12),b.h},"Fixed eta / pixel c or z0 / direct Q iteration",RGB(165,184,209),DT_CENTER);
        }
        auto v = state.cameras[state.plane];
        text(dc, {sc(12), cy(357 + panelExtra()), width - sc(24), sc(20)},
             spinMode() ? (spinAtlas() ? "k=["+num(state.spin.kMin,4)+", "+num(state.spin.kMax,4)+"] / tau=["+
                 num(state.spin.tauMin,4)+", "+num(state.spin.tauMax,4)+"] / drag pans; right-click selects"
                 : "Deterministic trajectories / burn-in >= averaging window / tangent growth is a finite-time estimate") :
             sphereMode() ? (state.generation==3 ? "Inverse-branch sample / post-burn-in points / conformal display / drag rotates cached states"
                                             : "Reduced-word sample / post-burn-in points / conformal display / drag rotates cached states") :
             quantumMode() ? "Finite-history density / 256 chains / post-burn-in samples / drag rotates cached states" :
             "Centre " + complexText({v.x, v.y}) + "    span " + num(v.span, 7) + "    " +
                 (state.family == 9 ? "x=log2(a), y=logit(r) / literal completed-pair Cayley diagnostic"
                  : (hilbertField() ? std::string("Direct Hilbert field    ") : std::to_string(state.iterations) + " iterations    ") +
                    (state.family ? (extraMode() ? frk_extra::extra_catalogue[state.generation-1].label : directQqrhMode() ? qqrhTitle(state.generation) : connectedLambdaMode() ? "BKQR lambda wrapper / finite RH host" : legacyLambdaMode() ? "Original BKQR lambda / finite ensemble diagnostics" : masterMode() ? masterTitle(state.generation) : tripleMode() ? momenttriple::title(state.generation) : quotientMode() ? "RH Jacobi quotient / direct rational iteration" : cubicMode() ? "RH moment cubic / direct iteration" : riemannMode() ? "RH finite-witness polynomial iteration" : projectMode() ? "Project complex iteration" : state.family == 11 ? "Spectral jet / direct rational iteration" : state.family == 10 ? (state.method == 1 ? "Halley + phase feedback" : "Newton + phase feedback") : (state.method == 1 ? "Halley garden" : "Newton garden"))
                                  : "R = c + alpha F/F(1)")),
             RGB(172, 183, 205));
        text(dc, {sc(12), cy(380 + panelExtra()), width - sc(24), sc(19)}, message, RGB(228, 205, 150));
        SelectObject(dc, smallFont);
        text(dc, {sc(12), cy(400 + panelExtra()), width - sc(24), sc(17)},
             directQqrhMode() ? "Right-click: select c / trace   J: Julia   C/B: seeded c portrait   K: Taylor N   P: palette   [ / ]: hue   , / .: bands" :
             connectedQqrhLambdaMode() ? "Square: drag / wheel / right-click readout   B: branch   C: outer H   K: Taylor N   P/M: metric   R: reset   Tab: QQRH channel" :
             connectedRhLambdaMode() ? "Square: drag / wheel / right-click readout   B: branch   C: outer H   K: Taylor N   P/M: metric   R: reset   Tab: RH channel" :
             legacyLambdaMode() ? "Square: drag to pan / wheel to zoom / right-click readout   F: host   B: branch   C: coordinates   P/M: metric   R: reset   Tab: catalogue" :
             masterMode() ? "T: preset bank   1-4: presets   V: controls   F: map   D: chart   I: coloring   J: Julia   C/B: c   A: lambda   K: degree" :
             hilbertField() ? "Drag: pan   Wheel: zoom   Right-click: inspect   I: inspector   P: palette   R: reset   H: formulas   Tab: catalogue" :
             extraMode() ? "J: starting z   C/B: parameter c   P: palette   R: reset   Right-click: select c / orbit   Wheel: zoom   Tab: catalogue" :
             directAtlasMode() ? "Drag: pan   Wheel: zoom   Right-click: select c / trace   J: Julia   C/B: atlas   1-4: presets   H: help   Tab: catalogue" :
             (projectMode() || riemannMode()) ? "Drag: pan   Wheel: zoom   Right-click: orbit   1-4: presets   S: BMP   H: help   Tab: catalogue   F11: fullscreen" :
             spinMode() ? (spinAtlas() ? "Drag: pan   Wheel: zoom   Right-click: select k/tau   E: selected sphere   B: atlas   J/C/A: views   H: help"
                 : state.generation==3 ? "Drag: rotate   Wheel: zoom   B: parameter atlas   J/C/A: views   N: new initial points   H: help   Tab: catalogue"
                                      : "Drag: rotate   Wheel: zoom   J/C/A: views   N: new initial points   E/B: samples   H: help   Tab: catalogue") :
             sphereMode() ? "Drag: rotate   Wheel: zoom   J/C/A: views   N: new sample   E/B: sample budget   S: BMP   H: help   Tab: catalogue" :
             quantumMode() ? "Drag: rotate   Wheel: zoom   J/C/A: views   N: new history   E/B: sample budget   S: BMP   H: help   Tab: catalogue" :
             state.family == 9
                 ? "Drag: synchronized pan   Wheel: synchronized zoom   Right-click: inspect a,r,C   S: BMP   H: help   Tab: renderers   F11: fullscreen"
                 : "Drag: pan   Wheel: zoom   Right-click: orbit   J/C/A: planes   E: explore   B: back   S: BMP   H: help   Tab: renderers   F11: fullscreen",
             RGB(127, 136, 160));
        SelectObject(dc, oldFont);
        presentPaint(target,dc);
    }
    int visibleEntries() const { return std::max(1, (height - sc(230)) / sc(64)); }
    Box catalogueRow(int i) const {
        return {sc(24), sc(120) + (i - catalogueScroll) * sc(64), width * 55 / 100 - sc(40), sc(57)};
    }
    Box previewBox() const {
        int x = width * 55 / 100 + sc(10), w = width - x - sc(26);
        return {x, sc(120), w, std::min(w * 3 / 4, std::max(sc(100), height - sc(485)))};
    }
    std::string entryDetail(const catalogue::Entry &e) const {
        if (!e.family)
            return "Original rational map / direct iteration";
        if(e.family==24)return std::string(qqrhTitle(e.generation))+" / direct iteration / seeded c portrait and selected Julia";
        if(e.family==23)return std::string(qqrhTitle(e.generation))+" / original BKQR outer lambda / independent inner eta";
        if(e.family==22)return std::string(masterTitle(e.generation))+" / original BKQR outer lambda / independent inner eta";
        if(e.family==21)return "Original BKQR lambda map / five hosts / three coordinate maps / nine sampled metrics";
        if(e.family==20) return std::string(masterTitle(e.generation))+" / GPU + CPU precision fallback";
        if(e.family==25) return std::string(frk_extra::extra_catalogue[e.generation-1].label)+
            " / native D3D11 GPU + parallel CPU fallback";
        if(e.family==19) return std::string(momenttriple::title(e.generation))+" / all-critical atlas and selected Julia plane";
        if(e.family==18) return "Five RH moments / Jacobi rational quotient / critical escape and Julia planes";
        if(e.family==17) return "Five RH moments / direct cubic / both-critical atlas and Julia planes";
        if(e.family==16) return std::string(riemannTitle(e.generation))+
            (e.generation>=4 ? " / fixed H2 moments / Newton + Halley" : " / finite witness / Newton + Halley");
        if(e.family==15) return std::string(projectTitle(e.generation))+" / 2D dynamics / CPU double";
        if (e.family == 14) return e.generation==3 ? "Twist/damping scan / finite-time stability / selected sphere exploration"
                                     : e.generation==2 ? "Analytic tangent dynamics / finite-time stretching / sphere and atlas"
                                                       : "Deterministic kicked spin / dissipative trajectories / cached sphere cloud";
        if (e.family == 13) return e.generation==3 ? "Degree-six rational map / cached inverse-branch sampling / sphere and atlas"
                                     : e.generation==2 ? "Separated circle reflections / Schottky limit set / cached reduced words"
                                                       : "Tangent circle reflections / Apollonian gasket / cached reduced words";
        if (e.family == 12) return "Born-weighted weak measurements / cached parallel CPU histories / sphere and equal-area views";
        if (e.family == 11)
            return e.generation == 2 ? "Three positive atoms / adjustable jets 1-3 / GPU + CPU deep zoom"
                                     : "Two positive spectral atoms / direct rational map / GPU + CPU deep zoom";
        auto d = garden::data(e.family, e.generation);
        if (e.family == 9)
            return "Literal completed-pair Cayley observable / 3 synchronized panels / CPU double";
        if (e.family == 10)
            return "Primitive chi mod 5 / Cayley sheet 3 / bounded phase-feedback Newton + Halley / CPU double";
        if (e.family == 7)
            return "Primitive chi mod 5 / Cayley sheet k=" + std::to_string(e.generation) + " / CPU double / Newton + Halley";
        if (e.family == 8)
            return "Primitive chi5 x chibar5 / completed pair / CPU double / Newton + Halley";
        if (e.family == 6)
            return "Analytic L jets / CPU double / Newton + Halley";
        if (e.family == 5)
            return "Cayley depth " + std::to_string(e.generation) + " / max degree " + std::to_string(d->degree) + " / Newton + Halley";
        if (e.family == 4)
            return std::to_string(d->transfer.states) + " states / spectral degree " +
                   std::to_string(d->degree) + " / Newton + Halley";
        if (e.family == 3)
            return "Depth " + std::to_string(e.generation) + " / " + std::to_string(d->degree) +
                   " valuation words / Newton + Halley";
        if (e.family == 1)
            return "Generation " + std::to_string(e.generation) + " / degree " + std::to_string(d->degree) +
                   " / Newton + Halley";
        return std::to_string(e.generation) + " x 1 x 1 cells / " + std::to_string(d->degree) +
               " bonds / Newton + Halley";
    }
    void paragraph(HDC dc, Box b, const std::string &value, COLORREF colour) {
        SetTextColor(dc, colour);
        SetBkMode(dc, TRANSPARENT);
        RECT r = b.rect();
        auto s = wide(value);
        DrawTextW(dc, s.c_str(), int(s.size()), &r, DT_LEFT | DT_WORDBREAK | DT_NOPREFIX);
    }
    void paintCatalogue(HDC dc) {
        fill(dc, {0, 0, width, height}, RGB(10, 11, 18));
        fill(dc, {0, 0, width, sc(94)}, RGB(20, 18, 30));
        auto old = SelectObject(dc, titleFont);
        text(dc, {sc(25), sc(20), width - sc(50), sc(32)}, "FRAKTALFAKTUR", RGB(239, 224, 255));
        SelectObject(dc, old);
        text(dc, {sc(26), sc(57), width - sc(52), sc(22)},
             "Choose a renderer. Every generation has its own entry.", RGB(171, 180, 204));
        std::string badge = std::to_string(entries.size()) + " RENDERERS  /  NATIVE WINDOWS";
        text(dc, {width - sc(360), sc(27), sc(334), sc(22)}, badge, RGB(185, 153, 241), DT_RIGHT);
        for (int i = catalogueScroll; i < int(entries.size()) && i < catalogueScroll + visibleEntries();
             i++) {
            auto b = catalogueRow(i);
            bool selected = i == selectedEntry;
            fill(dc, b, selected ? RGB(63, 44, 94) : RGB(24, 25, 36));
            border(dc, b, selected ? RGB(186, 155, 246) : RGB(56, 58, 76));
            text(dc, {b.x + sc(12), b.y + sc(8), sc(32), sc(34)}, (i < 9 ? "0" : "") + std::to_string(i + 1),
                 selected ? RGB(234, 218, 255) : RGB(133, 142, 166));
            text(dc, {b.x + sc(54), b.y + sc(6), b.w - sc(68), sc(22)}, entries[i].title, RGB(235, 231, 245));
            text(dc, {b.x + sc(54), b.y + sc(29), b.w - sc(68), sc(18)}, entryDetail(entries[i]),
                 selected ? RGB(198, 182, 224) : RGB(145, 155, 178));
        }
        auto box = previewBox();
        fill(dc, box, RGB(5, 6, 12));
        border(dc, box, RGB(104, 79, 135));
        if (frame && frame->generation == generation.load()) {
            BITMAPINFO bi{};
            bi.bmiHeader.biSize = 40;
            bi.bmiHeader.biWidth = frame->w;
            bi.bmiHeader.biHeight = -frame->h;
            bi.bmiHeader.biPlanes = 1;
            bi.bmiHeader.biBitCount = 32;
            StretchDIBits(dc, box.x + 1, box.y + 1, box.w - 2, box.h - 2, 0, 0, frame->w, frame->h,
                          frame->pixels.data(), &bi, DIB_RGB_COLORS, SRCCOPY);
        } else
            text(dc, box, "Rendering preview...", RGB(179, 159, 213), DT_CENTER);
        int y = box.y + box.h + sc(16);
        auto &e = entries[selectedEntry];
        text(dc, {box.x, y, box.w, sc(24)}, e.title, RGB(240, 226, 255));
        y += sc(33);
        paragraph(dc, {box.x, y, box.w, sc(62)}, e.description, RGB(178, 188, 209));
        y += sc(68);
        paragraph(dc, {box.x, y, box.w, sc(62)}, formula(state), RGB(198, 166, 241));
        y += sc(66);
        if(e.family==24)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                "Direct QQRH iteration with fixed complex eta. Scan c from seed zero, then select a parameter to explore its Julia plane and finite orbit traces.",RGB(163,171,194));
        else if(e.family==23)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                "Five QQRH hosts inside the original BKQR lambda wrapper. Inner eta and outer lambda are independent. Real positivity statements do not certify complex orbits.",RGB(163,171,194));
        else if(e.family==22)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                "The RH channel is the host inside the original BKQR lambda wrapper. Outer lambda scans the square; inner eta, Taylor degree and RH controls stay independent.",RGB(163,171,194));
        else if(e.family==21)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                "The original lambda-only ensemble map, with host, branch policy, coordinate map and all nine sampled metrics. Pan and zoom the square domain.",RGB(163,171,194));
        else if(e.family==25 && e.generation>=9)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                "Finite QQRH Hardy geometry: normalized coefficient vectors, shift correlations and angular derivative. N=12; radius R controls r=200R. H opens the formulas.",RGB(163,171,194));
        else if(e.family==25)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                "New mathematical dynamics, integrated with the original v1.8 GPU renderer and zoom controls. Finite trajectories are numerical pictures, not RH/Collatz proof.",RGB(163,171,194));
        else if(e.family==20)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                "Eighteen RH spectral maps. New entries 14-18: finite spectral echo, Laguerre determinant, dyadic holonomy, antipodal current and projective Schwarzian.",RGB(163,171,194));
        else if(e.family==19)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                "Three mechanisms from the moment data: conjugation, alternating phases, and pole feedback. Right-click a parameter to explore its Julia plane.",RGB(163,171,194));
        else if(e.family==18)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                "Iterate the quotient of two consecutive orthogonal polynomials. Scan both critical orbits, then right-click to explore the infinity basin boundary.",RGB(163,171,194));
        else if(e.family==17)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                "Scan both critical orbits of a cubic built from the finite RH Jacobi model. Right-click a parameter to open its Julia plane.",RGB(163,171,194));
        else if(e.family==16)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                e.generation>=4 ? "Explore Newton/Halley basins from five fixed RH-branch moments: orthogonal pencils, mirrored kernels, a finite Jacobi closure and Joukowski sheets."
                                : "Build a finite moment witness and explore its complex Newton/Halley basins. Compare direct nodes, the cubic pullback and reciprocal sheets.",RGB(163,171,194));
        else if(e.family==15)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                "One complex iteration per project. Explore finite root/cycle diagnostics, pan and zoom the boundary, and right-click to inspect an orbit.",RGB(163,171,194));
        else if (e.family == 14)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                      e.generation==3 ? "Scan twist and damping at fixed tilt. Select any parameter with a right-click and explore its trajectories on the sphere."
                                      : "Tilt, shear and damp a unit spin. Compare trajectory density with analytic tangent growth and local area change.",RGB(163,171,194));
        else if (e.family == 13)
            paragraph(dc,{box.x,y,box.w,sc(50)},
                      e.generation==3 ? "Explore backward orbits of z^3+lambda/z^3 on the Riemann sphere. Adjust complex lambda; compare the sphere, atlas and both hemispheres."
                                      : "Explore a reflection-group limit set directly on the sphere. Adjust the conformal display and inspect both hemispheres; Schottky adds circle separation.",RGB(163,171,194));
        else if (e.family == 12)
            paragraph(dc, {box.x,y,box.w,sc(50)},
                      "Explore conditioned measurement histories on a sphere. Compare their density with the decay of the ensemble mean; rotate without resampling.", RGB(163,171,194));
        else if (e.family == 11)
            paragraph(dc, {box.x,y,box.w,sc(50)},
                      e.generation == 2
                         ? "Three-atom spectral dynamics. Adjust the middle weight and choose the first, second or third jet to reshape the rational map."
                         : "Explore the rational dynamics of a spectral jet determinant. Change c, eta and atom separation d; compare with the eta=0 quadratic control.", RGB(163,171,194));
        else if (e.family == 9)
            paragraph(dc, {box.x, y, box.w, sc(50)},
                      "Literal completed-pair GRH closure observable specialized to primitive chi mod 5. Compare Theta phase, its logarithmic derivative, and C_chi(a,r) on one shared theorem-parameter plane.",
                      RGB(163, 171, 194));
        else if (e.family == 10)
            paragraph(dc, {box.x, y, box.w, sc(50)},
                      "Exact completed-L Cayley Newton/Halley jets plus a bounded vanishing residual-phase kick. Zeros stay fixed; only the surrounding basin dynamics are deformed.", RGB(163,171,194));
        else if (e.family == 2)
            paragraph(dc, {box.x, y, box.w, sc(50)},
                      "Finite graph connectivity; this is not an estimate of the "
                      "infinite-lattice threshold.",
                      RGB(163, 171, 194));
        else if (e.family == 8)
            paragraph(dc, {box.x, y, box.w, sc(50)},
                      "Target zero renders the completed chi5 x chibar5 pair. The two primitive zero families are combined into one mirror-coupled Newton/Halley garden.", RGB(163,171,194));
        else if (e.family == 7)
            paragraph(dc, {box.x, y, box.w, sc(50)},
                      "The critical line maps to Re(w)=0. Target zero renders the completed primitive L zero basins on one Cayley sheet; the vertical carrier is drawn in the canvas.", RGB(163,171,194));
        else if (e.family == 6)
            paragraph(dc, {box.x, y, box.w, sc(50)},
                      "Explore actual L-function zeros and complex level sets. The numerical window is bounded; these images are not RH/GRH certificates.", RGB(163,171,194));
        else if (e.family == 5)
            paragraph(dc, {box.x, y, box.w, sc(50)},
                      "The exact Cayley transmission kernel becomes a branching preimage tree. Move target or switch to Halley to reshape its basins.", RGB(163,171,194));
        else if (e.family == 4)
            paragraph(dc, {box.x, y, box.w, sc(50)},
                      "Complex feedback moves the eigenvalues of a finite restart operator. The fixed zero "
                      "factor is removed.",
                      RGB(163, 171, 194));
        else if (e.family == 3)
            paragraph(dc, {box.x, y, box.w, sc(50)},
                      "Move the complex embedding q to reshape the cylinder roots. "
                      "Colour follows the reached word.",
                      RGB(163, 171, 194));
        else if (e.family == 1)
            paragraph(dc, {box.x, y, box.w, sc(50)},
                      "Change complex branch weights a and b, or switch between "
                      "Newton and Halley.",
                      RGB(163, 171, 194));
        else
            paragraph(dc, {box.x, y, box.w, sc(50)},
                      "The original AMP renderer and its presets. Your version-1 "
                      ".view files still open.",
                      RGB(163, 171, 194));
        line(dc, sc(24), height - sc(102), width - sc(24), height - sc(102), RGB(59, 51, 77));
        button(dc, {sc(24), height - sc(83), sc(220), sc(36)}, "Open renderer [Enter]", true);
        button(dc, {sc(254), height - sc(83), sc(160), sc(36)}, "Load saved view");
        button(dc, {sc(424), height - sc(83), sc(90), sc(36)}, "Quit");
        text(dc, {sc(24), height - sc(37), width - sc(48), sc(20)},
             catalogueNotice ? message
                             : "Click to select / Enter or double-click to open / "
                               "Up and Down to browse / "
                               "Tab returns here",
             RGB(150, 159, 182));
        text(dc, {sc(530), height - sc(80), width - sc(556), sc(30)},
             hqBusy
                 ? "HQ export is rendering..."
                 : (frame && frame->generation == generation.load() ? frame->backend : "Preparing preview"),
             RGB(148, 162, 190), DT_RIGHT);
    }
    void selectEntry(int i) {
        i = std::clamp(i, 0, int(entries.size()) - 1);
        selectedEntry = i;
        catalogueNotice = false;
        if (i < catalogueScroll)
            catalogueScroll = i;
        if (i >= catalogueScroll + visibleEntries())
            catalogueScroll = i - visibleEntries() + 1;
        state = entryStates[i];
        orbitVisible = false;
        help = false;
        refresh();
        request();
    }
    void openEntry() {
        inCatalogue = false;
        activeEntry = selectedEntry;
        backValid = false;
        orbitVisible = false;
        help = false;
        frame.reset();
        refresh();
        request();
        SetWindowTextW(window, wide("FRAKTALFAKTUR / " + entries[activeEntry].title).c_str());
        notice("Tab / Renderers returns to the catalogue. Settings are remembered "
               "while this window is open.");
    }
    void showCatalogue() {
        if (inCatalogue)
            return;
        if (activeEntry >= 0 && activeEntry < int(entryStates.size()))
            entryStates[activeEntry] = state;
        inCatalogue = true;
        selectedEntry = std::max(0, activeEntry);
        help = false;
        orbitVisible = false;
        dragging = false;
        active = -1;
        ReleaseCapture();
        SetWindowTextW(window, title);
        selectEntry(selectedEntry);
    }
    void identifyState() {
        activeEntry = -1;
        for (size_t i = 0; i < entries.size(); i++)
            if (entries[i].family == state.family && entries[i].generation == state.generation) {
                activeEntry = int(i);
                entryStates[i] = state;
                selectedEntry = int(i);
                break;
            }
    }
    void catalogueClick(int x, int y, bool doubleClick = false) {
        for (int i = catalogueScroll; i < int(entries.size()) && i < catalogueScroll + visibleEntries(); i++)
            if (catalogueRow(i).has(x, y)) {
                selectEntry(i);
                if (doubleClick)
                    openEntry();
                return;
            }
        if (Box{sc(24), height - sc(83), sc(220), sc(36)}.has(x, y))
            openEntry();
        else if (Box{sc(254), height - sc(83), sc(160), sc(36)}.has(x, y))
            viewDialog(false);
        else if (Box{sc(424), height - sc(83), sc(90), sc(36)}.has(x, y))
            SendMessageW(window, WM_CLOSE, 0, 0);
    }
    void resetLabJuliaCamera() {
        state.cameras[1]=Camera{0,0,state.rhExperiment.chart==2 ? 12. : 4.};
    }
    C masterPhysicalPicked() const {
        return state.plane==1 ? rhmaster::chartForward(picked,state.rhExperiment.chart) : picked;
    }
    void cycleLab(int control) {
        auto &e=state.rhExperiment;
        if(control==0)e.map=(e.map+1)%rhmaster::MapCount;
        else if(control==1) {e.chart=(e.chart+1)%rhmaster::ChartCount;resetLabJuliaCamera();orbitVisible=false;}
        else e.coloring=(e.coloring+1)%rhmaster::ColouringCount;
        state.preset=-1;refresh();request();
        notice(std::string("Map: ")+rhmaster::mapName(e.map)+" / chart: "+rhmaster::chartName(e.chart)+
            " / color: "+rhmaster::colouringName(e.coloring)+". V toggles the lab sliders.");
    }
    void toggleLabControls() {
        labControls=!labControls;
        notice(labControls ? std::string("Lab bank: ")+focus_presets::bankName(focusBank)+". T switches banks; 1..4 select presets; V opens main controls."
                           : "Main controls. F map; D chart; I coloring; T bank; V returns to lab sliders.");
        InvalidateRect(window,nullptr,FALSE);
    }
    void cycleFocusBank() {
        focusBank=(focusBank+1)%focus_presets::BankCount;
        // Selecting a bank only changes its available presets; retain the current image.
        state.preset=-1;
        notice(std::string("Preset bank: ")+focus_presets::bankName(focusBank)+". Press 1..4 in Lab controls to explore; T changes bank.");
        InvalidateRect(window,nullptr,FALSE);
    }
    void applyLabPreset(int index) {
        focus_presets::apply(state,focusBank,index);
        orbitVisible=false;refresh();request();
        notice(std::string(focus_presets::bankName(focusBank))+": "+focus_presets::presetName(focusBank,index)+
            " / "+rhmaster::colouringName(state.rhExperiment.coloring)+". T changes preset bank.");
    }
    bool labPresetSelected(int index) const {
        if(state.preset!=index)return false;
        auto candidate=state;
        focus_presets::apply(candidate,focusBank,index);
        return candidate.c==state.c && candidate.alpha==state.alpha &&
            rhmaster::sameExperiment(candidate.rhExperiment,state.rhExperiment);
    }
    void changePlane(int p) {
        if(legacyLambdaMode() && p!=2)return;
        if(extraMode() && (p<0 || p>1))return;
        if(extraMode() && p==0 && frk_extra::seed_at_pixel(state.generation+17)) {
            notice("This root-map or entire-map field has no variable c-plane.");return;
        }
        if(masterMode() && (p<0 || p>2)) return;
        if(directAtlasMode() && !masterMode() && (p<0 || p>1)) return;
        if((projectMode() || riemannMode()) && p!=0) return;
        if (p == state.plane)
            return;
        state.plane = p;
        orbitVisible = false;
        refresh();
        request();
    }
    Box spinPlotBox() const {
        int w=width,h=canvasH;
        if (frame && frame->state.family==14 && frame->state.generation==3 && frame->state.plane==0) {
            w=frame->w;h=frame->h;
        }
        const auto b=spin::AtlasRenderer::plotRect(w,h);
        const int left=int(std::lround(double(b.x)*width/std::max(1,w)));
        const int top=int(std::lround(double(b.y)*canvasH/std::max(1,h)));
        const int right=int(std::lround(double(b.x+b.width)*width/std::max(1,w)));
        const int bottom=int(std::lround(double(b.y+b.height)*canvasH/std::max(1,h)));
        return {left,top,std::max(1,right-left),std::max(1,bottom-top)};
    }
    static void shiftRange(double &lo, double &hi, double delta, double maximum) {
        const double span=std::clamp(hi-lo,1.000001e-6,maximum);
        lo=std::clamp(lo+delta,0.,maximum-span);
        hi=lo+span;
    }
    void panSpinAtlas(double dx, double dy) {
        shiftRange(state.spin.kMin,state.spin.kMax,dx,20);
        shiftRange(state.spin.tauMin,state.spin.tauMax,dy,2);
        request();
    }
    void resetCamera() {
        if(masterMode() && state.plane==1) {resetLabJuliaCamera();request();return;}
        if (spinAtlas()) {
            const auto t=defaults(14,3).spin;
            state.spin.kMin=t.kMin;state.spin.kMax=t.kMax;
            state.spin.tauMin=t.tauMin;state.spin.tauMax=t.tauMax;
            request();return;
        }
        if (sampleMode()) {
            state.quantum.yaw=.4; state.quantum.pitch=.2; state.quantum.zoom=1;
            request(); return;
        }
        state.cameras[state.plane] = defaults(state.family, state.generation).cameras[state.plane];
        request();
    }
    void explore() {
        if(legacyLambdaMode())return;
        if(extraMode()) {changePlane(1);notice("Starting-state plane / selected c="+complexText(state.c));return;}
        if(directQqrhMode()) {changePlane(1);notice("Selected c="+complexText(state.c)+", eta="+complexText(state.alpha)+". C / B returns to the seeded c portrait.");return;}
        if(masterMode()) {changePlane(1);notice("Selected c="+complexText(state.c)+", lambda="+complexText(state.alpha)+". C / B: c portrait; A: lambda plane.");return;}
        if(directAtlasMode()) {changePlane(1);notice(std::string(tripleMode() ? "Selected moment c=" : quotientMode() ? "Selected quotient c=" : "Selected cubic c=")+complexText(state.c)+". C / B returns to the atlas.");return;}
        if (spinMode() && state.generation==3) {
            changePlane(1);
            notice("Selected spin parameters: k="+num(state.spin.twist,5)+", tau="+num(state.spin.damping,5)+". B returns to the atlas.");
            return;
        }
        if (sampleMode()) { setValue(2,state.quantum.samplePower+1); return; }
        if (state.plane == 0)
            return;
        back = state;
        backValid = true;
        state.plane = 0;
        orbitVisible = false;
        refresh();
        request();
        notice("Exploring the selected parameter. B / Back restores the atlas.");
    }
    void goBack() {
        if(legacyLambdaMode())return;
        if(directAtlasMode()) {changePlane(0);return;}
        if (spinMode() && state.generation==3) {changePlane(0);return;}
        if (sampleMode()) { setValue(2,state.quantum.samplePower-1); return; }
        if (!backValid)
            return;
        state = back;
        backValid = false;
        orbitVisible = false;
        refresh();
        request();
    }
    fs::path output(const wchar_t *ext) {
        fs::create_directories(root / L"renders");
        return root / L"renders" /
               (wide("fraktal_" + std::to_string(state.family) + "_g" + std::to_string(state.generation) +
                     "_" + stamp()) +
                ext);
    }
    void saveBmp() {
        if (!frame || frame->generation != generation.load() || frame->stride != 1 ||
            (masterMode() && !samePaintScene(frame->state,state))) {
            notice("Wait for FULL RES before saving, or use Save HQ.");
            return;
        }
        auto p = output(L".bmp");
        writeBmp(p, frame->pixels, frame->w, frame->h);
        State saved = frame->state;
        if (saved.family == 11) {
            // Seed and guide toggles do not alter pixels in the x-start plane.
            saved.seed = state.seed;
            saved.method = state.method;
        }
        metadata(p, saved, frame->w, frame->h, 1, frame->backend, &frame->quantumStats, &frame->sphereStats, &frame->spinStats, frame->projectModel.get(), frame->riemannModel.get(), frame->cubicModel.get(), frame->quotientModel.get(), frame->tripleModel.get(), frame->masterModel.get(), frame->qqrhDirectModel.get());
        notice("BMP saved: " + utf8(p.filename().wstring()));
    }
    void saveHQ() {
        if (hqBusy)
            return;
        Job j;
        j.generation = ++exportGeneration;
        j.state = state;
        j.w = std::max(1920, width);
        j.h = std::max(1, int(std::round(double(canvasH) * j.w / width)));
        j.output = output(L"_HQ.bmp");
        {
            std::lock_guard<std::mutex> lock(mutex);
            exportJob = j;
            hqPending = true;
        }
        hqBusy = true;
        cv.notify_one();
        notice("Rendering HQ at " + std::to_string(j.w) + " x " + std::to_string(j.h) +
               (legacyLambdaMode() ? ", larger canvas with the same grid/seed settings. Esc cancels." : spinAtlas() ? ", denser parameter grid. Esc cancels." : sampleMode() ? ", up to 4x retained samples. Esc cancels." : ", 2 x 2 samples. Esc cancels."));
    }
    void newSampleHistory() {
        state.quantum.randomSeed = state.quantum.randomSeed==0xffffffffu ? 1u : state.quantum.randomSeed+1u;
        state.preset=-1;
        request();
        notice(std::string(spinMode() ? "New initial points; deterministic map unchanged. Seed = " : sphereMode() ? "New spherical sample; replay seed = " : "New measurement history; replay seed = ")+std::to_string(state.quantum.randomSeed));
    }
    void saveOrbit() {
        if(legacyLambdaMode()) {notice("This lambda-only map exports ensemble settings with each image; no single orbit is selected.");return;}
        if (sampleMode()) {
            if (!frame || frame->generation!=generation.load() || frame->stride!=1) {
                notice("Wait for the completed sample before saving statistics."); return;
            }
            const auto p=output(L"_stats");
            metadata(p,frame->state,frame->w,frame->h,1,frame->backend,&frame->quantumStats,&frame->sphereStats,&frame->spinStats);
            notice("Statistics saved: "+utf8(p.filename().wstring())+".txt");
            return;
        }
        if (!orbitVisible) {
            notice("Right-click the canvas to trace an orbit first.");
            return;
        }
        auto p = output(L"_orbit.csv");
        std::ofstream f(p);
        if (!f)
            throw std::runtime_error("Cannot write orbit CSV");
        f << "orbit,step,z_real,z_imag,is_infinity,X_real,X_imag,Y_real,Y_imag\n" << std::setprecision(17);
        auto save = [&](const std::vector<Q> &trail, const char *name) {
            for (size_t i = 0; i < trail.size(); i++) {
                Q q = trail[i];
                C z;
                bool finite = affine(q, z);
                f << name << ',' << i << ',';
                if (finite)
                    f << z.real();
                f << ',';
                if (finite)
                    f << z.imag();
                f << ',' << (finite ? 0 : 1) << ',' << q.x.real() << ',' << q.x.imag() << ',' << q.y.real()
                  << ',' << q.y.imag() << '\n';
            }
        };
        save(trailA, (directAtlasMode() || state.plane == 0)         ? "picked_z"
                     : (state.family || state.orbitMode == 3) ? "custom_seed"
                                                              : "critical_plus");
        save(trailB, "critical_minus");
        f.close();
        if (!f)
            throw std::runtime_error("Orbit write failed");
        writeView(p.wstring() + L".view", state);
        std::ofstream info(fs::path(p.wstring() + L".txt"));
        info << std::setprecision(17) << "FRAKTALFAKTUR orbit: " << formula(state)
             << "\nFamily: " << familyName(state.family) << "\nGeneration: " << state.generation
             << "\nc = " << state.c;
        if(!directQqrhMode())info << "\nalpha = " << state.alpha << "\nSpectral d = " << state.spectralD
            << "\nMiddle weight = " << state.spectralWeight << "\nJet order = " << state.spectralOrder;
        info << "\nPicked z = " << picked << (directQqrhMode()
                 ? "\nFinite direct-Q orbit prefix; iteration limits and numerical limits remain unresolved.\n"
                 : masterMode()
                 ? "\nFinite rational-feedback orbit prefix; iteration limits and numerical limits remain unresolved.\n"
                 : tripleMode()
                 ? "\nFinite map-orbit prefix; a return step contains two quadratic phases in generation 2. Homogeneous endpoint recorded separately.\n"
                 : quotientMode()
                 ? "\nFinite direct-rational orbit prefix; stops before infinity. The homogeneous endpoint is recorded separately.\n"
                 : cubicMode()
                 ? "\nFinite direct-cubic orbit prefix; budget-limited orbits remain unresolved.\n"
                 : (projectMode() || riemannMode())
                 ? "\nFinite orbit prefix only; poles and numerical cutoffs end the trace. No infinite tail or escape is inferred.\n"
                 : state.family == 11
                 ? "\nFinite orbit prefix only; infinity after a pole/escape is not appended.\n"
                 : "\nInfinity is explicitly flagged; homogeneous coordinates preserve the point.\n");
        if(directQqrhMode() && diagnostic.qqrhDirectModel) {
            const auto result=diagnostic.qqrhDirectModel->orbit(picked,state.c,state.iterations);
            info << "Plane: initial state z / selected c / fixed complex eta\n";
            directQqrhModelMetadata(info,diagnostic.qqrhDirectModel.get(),state.c);
            info << "Orbit status: " << qqrhdirect::statusName(result.status)
                 << "\nIteration budget: " << state.iterations
                 << "\nExecuted steps: " << result.iterations << "\nEndpoint: " << result.endpoint
                 << "\nDetected finite-tail near-period: " << result.detectedPeriod
                 << "\nAdjacent block comparisons: " << result.periodRepeats
                 << "\nRelative period residual: " << result.periodResidual
                 << "\nDetection tests periods up to 16 using four full adjacent-block comparisons across five stored blocks, at relative residual <=1e-7."
                 << "\nPeriod zero means no near-period detected; a detected near-period is not a certified stable cycle.\n";
        }
        if(riemannMode()) {
            info << "Moment order: H" << state.riemann.order
                 << (state.generation>=4 ? "\nH2 anchor (witness diagnostics only): " : "\nH2 anchor: ") << state.riemann.anchor
                 << "\nMethod: " << (state.method ? "Halley" : "Newton") << '\n';
            riemannModelMetadata(info,diagnostic.riemannModel.get());
        }
        if(cubicMode() && diagnostic.cubicModel) {
            const auto result=diagnostic.cubicModel->orbit(picked,state.c,state.iterations);
            info << "Plane: initial state z / selected c\nIteration: z <- z^3-A*z+c\nA: " << diagnostic.cubicModel->A << "\nCritical points: +/- " << diagnostic.cubicModel->critical
                 << "\nSelected-c escape radius: " << diagnostic.cubicModel->escapeRadius(state.c)
                 << "\nOrbit status: " << momentcubic::statusName(result.status)
                 << "\nExecuted steps: " << result.iterations << "\nEndpoint: " << result.endpoint << '\n';
            info << "BEGIN FIXED MOMENT SOURCE (not the iterated map)\n";
            riemannModelMetadata(info,diagnostic.cubicModel->source.get());
            info << "END FIXED MOMENT SOURCE\n";
        }
        if(masterMode() && diagnostic.masterModel) {
            const auto result=diagnostic.masterModel->orbitAtLambda(masterPhysicalPicked(),state.c,state.alpha,state.iterations);
            info << "Plane: Julia chart / selected c and complex lambda\n"
                 << "Picked chart coordinate: " << picked << "\nPicked physical state: " << masterPhysicalPicked()
                 << "\nCSV z and homogeneous coordinates record the displayed chart trace; the reported orbit endpoint is physical.\n";
            masterModelMetadata(info,diagnostic.masterModel.get(),state.c);
            info << "Orbit status: " << rhmaster::statusName(result.status)
                 << "\nExecuted steps: " << result.iterations << "\nEndpoint: " << result.endpoint
                 << "\nDetected finite-tail near-period: " << result.detectedPeriod
                 << "\nAdjacent block comparisons: " << result.periodRepeats
                 << "\nRelative period residual: " << result.periodResidual
                 << "\nDetection tests periods up to 16 using four full adjacent-block comparisons across five stored blocks, at relative residual <=1e-7."
                 << "\nPeriod zero means no near-period detected; a detected near-period is not a certified stable cycle.\n";
        }
        if(tripleMode() && diagnostic.tripleModel) {
            const auto result=diagnostic.tripleModel->orbit(picked,state.c,state.iterations);
            info << "Plane: initial state z / selected c\n";
            tripleModelMetadata(info,diagnostic.tripleModel.get(),state.c);
            info << "Orbit status: " << momenttriple::statusName(result.status)
                 << "\nExecuted steps: " << result.iterations
                 << "\nEndpoint is infinity: " << (result.atInfinity ? "yes" : "no")
                 << "\nEndpoint affine representable: " << (result.endpointRepresentable ? "yes" : "no")
                 << "\nHomogeneous endpoint [Z:W]: [" << result.endpointZ << ':' << result.endpointW << "]\n";
            if(result.atInfinity) info << "Endpoint: infinity\n";
            else if(result.endpointRepresentable) info << "Endpoint: " << result.endpoint << '\n';
            else info << "Endpoint: finite projective point; affine double is not representable\n";
            info << "Finite nonescape within the iteration budget remains unresolved.\n";
        }
        if(quotientMode() && diagnostic.quotientModel) {
            const auto result=diagnostic.quotientModel->orbit(picked,state.c,state.iterations);
            info << "Plane: initial state z / selected c\n";
            quotientModelMetadata(info,diagnostic.quotientModel.get(),state.c);
            info << "Orbit status: " << momentquotient::statusName(result.status)
                 << "\nExecuted steps: " << result.iterations
                 << "\nEndpoint is infinity: " << (result.atInfinity ? "yes" : "no")
                 << "\nEndpoint affine representable: " << (result.endpointRepresentable ? "yes" : "no")
                 << "\nHomogeneous endpoint [Z:W]: [" << result.endpointZ << ':' << result.endpointW << "]\n";
            if(result.atInfinity) info << "Endpoint: infinity\n";
            else if(result.endpointRepresentable) info << "Endpoint: " << result.endpoint << '\n';
            else info << "Endpoint: finite projective point; affine double is not representable\n";
            info << "Finite nonescape within the iteration budget remains unresolved.\n";
        }
        if(!info) throw std::runtime_error("Orbit metadata write failed");
        notice("Orbit saved: " + utf8(p.filename().wstring()));
    }
    void viewDialog(bool save) {
        fs::create_directories(root / L"views");
        wchar_t path[32768]{};
        auto initial = (root / L"views").wstring();
        if (save)
            lstrcpyW(path, L"my_fraktal_view.view");
        OPENFILENAMEW o{};
        o.lStructSize = sizeof(o);
        o.hwndOwner = window;
        o.lpstrFilter = L"Fraktalfaktur view (*.view)\0*.view\0\0";
        o.lpstrFile = path;
        o.nMaxFile = 32768;
        o.lpstrInitialDir = initial.c_str();
        o.lpstrDefExt = L"view";
        o.Flags = OFN_NOCHANGEDIR | OFN_PATHMUSTEXIST | (save ? OFN_OVERWRITEPROMPT : OFN_FILEMUSTEXIST);
        o.lpstrTitle = save ? L"Save renderer view" : L"Load renderer view";
        if (save ? GetSaveFileNameW(&o) : GetOpenFileNameW(&o)) {
            if (save) {
                writeView(path, state);
                notice("View saved.");
            } else {
                State loaded = readView(path);
                if (!inCatalogue && activeEntry >= 0)
                    entryStates[activeEntry] = state;
                state = loaded;
                identifyState();
                inCatalogue = false;
                frame.reset();
                SetWindowTextW(window, wide("FRAKTALFAKTUR / " + familyName(state.family)).c_str());
                backValid = false;
                orbitVisible = false;
                refresh();
                request();
                notice("View loaded: " + utf8(fs::path(path).filename().wstring()));
            }
        }
    }
    void toggleFullscreen() {
        if (!full) {
            oldStyle = GetWindowLongW(window, GWL_STYLE);
            placement.length = sizeof(placement);
            GetWindowPlacement(window, &placement);
            MONITORINFO m{};
            m.cbSize = sizeof(m);
            GetMonitorInfoW(MonitorFromWindow(window, MONITOR_DEFAULTTONEAREST), &m);
            SetWindowLongW(window, GWL_STYLE, oldStyle & ~WS_OVERLAPPEDWINDOW);
            SetWindowPos(window, HWND_TOP, m.rcMonitor.left, m.rcMonitor.top,
                         m.rcMonitor.right - m.rcMonitor.left, m.rcMonitor.bottom - m.rcMonitor.top,
                         SWP_FRAMECHANGED);
            full = true;
        } else {
            SetWindowLongW(window, GWL_STYLE, oldStyle);
            SetWindowPlacement(window, &placement);
            SetWindowPos(window, nullptr, 0, 0, 0, 0,
                         SWP_NOMOVE | SWP_NOSIZE | SWP_NOZORDER | SWP_FRAMECHANGED);
            full = false;
        }
    }
    void click(int x, int y) {
        if (inCatalogue) {
            catalogueClick(x, y);
            return;
        }
        if (Box{width - sc(165), cy(4), sc(153), sc(23)}.has(x, y)) {
            showCatalogue();
            return;
        }
        if (y < canvasH) {
            if (help) {
                help = false;
                InvalidateRect(window, nullptr, FALSE);
                return;
            }
            if(legacyLambdaMode() && !legacyLambdaBox().has(x,y))return;
            if (spinAtlas()) {
                dragSpinPlot=spinPlotBox();
                if (!dragSpinPlot.has(x,y)) return;
                dragSpin=state.spin;
            }
            dragging = true;
            dragX = x;
            dragY = y;
            dragCamera = sampleMode() ? Camera{state.quantum.yaw,state.quantum.pitch,state.quantum.zoom} : state.cameras[state.plane];
            SetCapture(window);
            return;
        }
        for (int i = 0; i < 4; i++)
            if (row(i, 4, 49).has(x, y)) {
                if(masterMode() && labControls) {applyLabPreset(i);return;}
                preset(state, i);
                orbitVisible = false;
                refresh();
                request();
                return;
            }
        if (!state.family) {
            for (int i = 0; i < 7; i++)
                if (row(i, 7, 81).has(x, y)) {
                    if (i < 3)
                        changePlane(i);
                    else if (state.plane != 0) {
                        state.orbitMode = i - 3;
                        refresh();
                        request();
                    }
                    return;
                }
        } else {
            for (int i = 0; i < 6; i++)
                if (row(i, 6, 81).has(x, y)) {
                    if(directQqrhMode()) {
                        if(i<2)changePlane(i);
                        else if(i==2) {state.method=(state.method+1)%3;state.preset=-1;refresh();request();}
                        else if(i==5) {help=!help;InvalidateRect(window,nullptr,FALSE);}
                    } else if(connectedLambdaMode()) {
                        if(i==2)state.orbitMode=(state.orbitMode+1)%3;
                        else if(i==3)state.spectralOrder=state.spectralOrder%3+1;
                        else if(i==4)state.method=(state.method+1)%3;
                        else if(i==5) {help=!help;InvalidateRect(window,nullptr,FALSE);}
                        if(i>=2 && i<=4) {state.preset=-1;refresh();request();}
                    } else if(legacyLambdaMode()) {
                        if(i==1)setValue(4,(state.method+1)%5);
                        else if(i==2)setValue(5,(state.orbitMode+1)%3);
                        else if(i==3)setValue(6,state.spectralOrder%3);
                        else if(i==4) {state.palette=(state.palette+1)%9;request();}
                        else if(i==5) {help=!help;InvalidateRect(window,nullptr,FALSE);}
                    } else if(masterMode()) {
                        if(i==4)toggleLabControls();
                        else if(labControls && i<3)cycleLab(i);
                        else if(i<3)changePlane(i);
                        else if(i==3 && labControls)cycleFocusBank();
                        else if(i==3) {state.spectralOrder=state.spectralOrder%3+1;state.preset=-1;refresh();request();}
                        else if(i==5) {help=!help;InvalidateRect(window,nullptr,FALSE);}
                    } else if(directAtlasMode()) {
                        if(i<2) changePlane(i);
                        else if(i==5) {help=!help;InvalidateRect(window,nullptr,FALSE);}
                    } else if(riemannMode()) {
                        if(i==0) changePlane(0);
                        else if(i==3 || i==4) {state.method=i-3;state.preset=-1;refresh();request();}
                        else if(i==5) {help=!help;InvalidateRect(window,nullptr,FALSE);}
                    } else if(projectMode()) {
                        if(i==0) changePlane(0);
                        else if(i==5) {help=!help;InvalidateRect(window,nullptr,FALSE);}
                    } else if (sampleMode()) {
                        if(i<3) changePlane(i);
                        else if(i==3) {state.method=1-state.method;request();}
                        else if(i==4) newSampleHistory();
                        else {help=!help;InvalidateRect(window,nullptr,FALSE);}
                    } else if (state.family == 9) {
                        if (i == 5) { help = !help; InvalidateRect(window, nullptr, FALSE); }
                    } else if (state.family == 11) {
                        if (i < 3) changePlane(i);
                        else if (i == 3) { state.method = 1-state.method; InvalidateRect(window,nullptr,FALSE); }
                        else if (i == 4) { state.seed = C(0,0); refresh(); request(); notice("Parameter seed reset to zero."); }
                        else { help = !help; InvalidateRect(window,nullptr,FALSE); }
                    } else if (i < 3)
                        changePlane(i);
                    else if (i < 5) {
                        state.method = i - 3;
                        refresh();
                        request();
                    } else {
                        help = !help;
                        InvalidateRect(window, nullptr, FALSE);
                    }
                    return;
                }
        }
        for (int i = 0; i < 3; i++)
            if (row(i, 3, 113).has(x, y)) {
                state.palette = legacyLambdaMode() ? (state.palette+(i==0 ? 8 : 1))%9 : i;
                request(true);
                return;
            }
        for (int i = 0; i < 8; i++)
            if (row(i, 8, 146, 29).has(x, y)) {
                switch (i) {
                case 0:
                    resetCamera();
                    break;
                case 1:
                    state = defaults(state.family, state.generation);
                    orbitVisible = false;
                    backValid = false;
                    refresh();
                    request();
                    break;
                case 2:
                    saveBmp();
                    break;
                case 3: {
                    fs::create_directories(root / L"renders");
                    auto result = ShellExecuteW(window, L"open", (root / L"renders").c_str(), nullptr,
                                                nullptr, SW_SHOWNORMAL);
                    if (INT_PTR(result) <= 32)
                        notice("Could not open renders folder.");
                    break;
                }
                case 4:
                    saveOrbit();
                    break;
                case 5:
                    if (sampleMode()) {newSampleHistory();break;}
                    orbitVisible = false;
                    trailA.clear();
                    trailB.clear();
                    InvalidateRect(window, nullptr, FALSE);
                    break;
                case 6:
                    help = !help;
                    InvalidateRect(window, nullptr, FALSE);
                    break;
                case 7:
                    SendMessageW(window, WM_CLOSE, 0, 0);
                    break;
                }
                return;
            }
        for (int i = 0; i < 5; i++)
            if (row(i, 5, 182, 29).has(x, y)) {
                switch (i) {
                case 0:
                    explore();
                    break;
                case 1:
                    goBack();
                    break;
                case 2:
                    viewDialog(true);
                    break;
                case 3:
                    viewDialog(false);
                    break;
                case 4:
                    saveHQ();
                    break;
                }
                return;
            }
        for (auto s : sliders())
            if (s.box.has(x, y) && enabled(s.id)) {
                active = s.id;
                SetCapture(window);
                moveSlider(x);
                return;
            }
    }
    void zoom(int x, int y, double factor) {
        if(legacyLambdaMode()) {
            const auto b=legacyLambdaBox();
            if(!b.has(x,y))return;
            auto &v=state.cameras[2];
            const C before=legacyLambdaWorld(x,y);
            v.span=std::clamp(v.span*factor,1e-10,1e5);
            const C after=legacyLambdaWorld(x,y);
            v.x+=before.real()-after.real();v.y+=before.imag()-after.imag();
            request();return;
        }
        if (spinAtlas()) {
            const auto b=spinPlotBox();
            if (!b.has(x,y)) return;
            const double fx=std::clamp(double(x-b.x)/b.w,0.,1.);
            const double fy=std::clamp(double(y-b.y)/b.h,0.,1.);
            auto &t=state.spin;
            const double k=t.kMin+fx*(t.kMax-t.kMin),tau=t.tauMax-fy*(t.tauMax-t.tauMin);
            const double ks=std::clamp((t.kMax-t.kMin)*factor,1.000001e-6,20.);
            const double ts=std::clamp((t.tauMax-t.tauMin)*factor,1.000001e-6,2.);
            t.kMin=std::clamp(k-fx*ks,0.,20-ks);t.kMax=t.kMin+ks;
            t.tauMin=std::clamp(tau-(1-fy)*ts,0.,2-ts);t.tauMax=t.tauMin+ts;
            request();return;
        }
        if (sampleMode()) {state.quantum.zoom=std::clamp(state.quantum.zoom/factor,.4,8.0);request();return;}
        auto &v = state.cameras[state.plane];
        int mapW = state.family == 9 ? std::max(1, width / 3) : width;
        int mapX = state.family == 9 ? (x % mapW) : x;
        C before = world(v, mapX, y, mapW, canvasH);
        v.span = std::clamp(v.span * factor, 1e-10, 1e5);
        C after = world(v, mapX, y, mapW, canvasH);
        v.x += before.real() - after.real();
        v.y += before.imag() - after.imag();
        if (state.family == 9) state.cameras[1] = state.cameras[2] = v;
        request();
    }
    void wheel(int x, int y, int delta) {
        if (inCatalogue) {
            selectEntry(selectedEntry - (delta > 0 ? 1 : -1));
            return;
        }
        if (y < canvasH) {
            zoom(x, y, std::pow(.8, delta / double(WHEEL_DELTA)));
            return;
        }
        for (auto s : sliders())
            if (s.box.has(x, y) && enabled(s.id)) {
                double step = directQqrhMode() ? (s.id==2 ? 16 : s.id==8 ? 1 : s.id==9 ? (state.generation==6 ? .025 : .1) : s.id==6 ? .005 : (s.id==3 || s.id==4) ? .025 : .01) :
                    connectedLambdaMode() ? (s.id==2 ? 16 : s.id==11 ? .1 : s.id>=8 ? 1 : s.id==6 ? .005 : (s.id==3 || s.id==4) ? .025 : .01) :
                    legacyLambdaMode() ? (s.id==2 ? 16 : s.id>=3 ? 1 : .01) :
                    masterMode() && labControls ? (s.id==0 ? 1 : s.id==2 ? 16 : s.id==6 ? .05 : s.id==8 ? .1 : .025) :
                    masterMode() ? (s.id==2 ? 16 : s.id==8 ? 1 : s.id==6 ? .005 : (s.id==3 || s.id==4) ? .025 : .01) :
                    riemannMode() ? (s.id==0 ? 1 : s.id==1 ? .02 : s.id==2 ? 16 :
                    s.id==5 ? .005 : s.id==8 ? .025 : .01) : projectMode() ? (s.id==2 ? 16 : s.id==5 ? .005 : s.id==8 ? .025 :
                    state.generation==1 && s.id==6 ? 2 : state.generation==1 && s.id==7 ? .25 :
                    state.generation==2 || (state.generation==3 && s.id==0) ? 1 :
                    state.generation==3 && s.id==1 ? .00001 : state.generation==5 ? .001 : .01) : spinMode() ? (s.id==2 || s.id==4 || s.id==10 || s.id==11 ? 1 :
                    s.id==9 ? .1 : s.id==8 ? .05 : s.id==3 ? .005 : .01) : sphereMode() ? (s.id==2 || s.id==11 ? 1 : s.id==8 ? .1 : s.id==7 ? .05 : s.id==0 ? .001 : .01) :
                    quantumMode() ? (s.id==2 || s.id==10 || s.id==11 ? 1 : s.id==9 ? .1 : s.id==8 ? .05 : .01) :
                    s.id == 10 ? 1 : s.id == 2 ? 16 : s.id == 5 ? .005 : s.id == 8 ? .025 : .01;
                setValue(s.id, value(s.id) + delta / double(WHEEL_DELTA) * step);
                return;
            }
    }
    void rightClick(int x, int y) {
        if (inCatalogue)
            return;
        if (y >= canvasH || help)
            return;
        if(legacyLambdaMode()) {
            if(!legacyLambdaBox().has(x,y))return;
            const C lambda=legacyLambdaWorld(x,y);
            notice(std::string(connectedLambdaMode() ? "outer lambda=" : "lambda=")+complexText(lambda)+
                (connectedLambdaMode() ? " / inner eta="+complexText(state.alpha) : std::string{})+" / "+bkqrlambda::branchName(state.orbitMode)+
                " / "+bkqrlambda::couplingName(state.spectralOrder-1)+" / "+legacyMetricTitle(state.palette)+". The lambda map remains active.");
            return;
        }
        if(extraMode()) {
            if(!frame || frame->generation!=generation.load() || frame->state.family!=25 ||
               frame->state.generation!=state.generation || frame->state.plane!=state.plane) {
                notice("Wait for the new fractal preview to finish before selecting an orbit.");return;
            }
            const C clicked=world(state.cameras[state.plane],x,y,width,canvasH);
            if(hilbertField()) {
                hilbertPicked=clicked;hilbertInspectorVisible=true;picked=clicked;orbitVisible=false;refresh();
                const auto& model=qqrh_hilbert::cached(state.alpha.real());
                notice("z="+complexText(clicked)+" / overlap="+num(model.coherence(clicked),7)+
                    " / angle="+num(model.projective_distance(clicked),5)+" rad / g="+
                    num(model.projective_metric(clicked),6));return;
            }
            if(state.plane==0 && !frk_extra::seed_at_pixel(state.generation+17)) {
                state.c=clicked;state.preset=-1;changePlane(1);
                notice("Selected c="+complexText(state.c)+" / opening Julia basin.");
            } else {
                picked=clicked;orbitVisible=true;refresh();
                frk_extra::Params q;q.lambda=state.alpha.real();q.julia_c=state.c;
                q.julia=true;q.iterations=state.iterations;
                const auto orbit=frk_extra::evaluate(state.generation+17,picked,q);
                notice("Orbit start="+complexText(picked)+" / steps="+std::to_string(orbit.iterations)+
                       " / "+(orbit.converged?"root reached":orbit.escaped?"escaped":"nonescape unresolved")+
                       " / endpoint="+complexText(orbit.final_z));
            }
            return;
        }
        if(directQqrhMode()) {
            if(!diagnostic.qqrhDirectModel || !diagnostic.qqrhDirectModel->validModel) {
                notice("Model unavailable: "+(diagnostic.qqrhDirectModel ? diagnostic.qqrhDirectModel->statusText : std::string("missing model")));
                return;
            }
            if(!frame || frame->generation!=generation.load() || frame->state.family!=24 || frame->state.generation!=state.generation || frame->state.plane!=state.plane) {
                notice("Wait for the current direct QQRH preview before selecting a parameter or orbit.");return;
            }
            const C clicked=world(state.cameras[state.plane],x,y,width,canvasH);
            if(state.plane==0) {
                if(!std::isfinite(clicked.real()) || !std::isfinite(clicked.imag()) || std::abs(clicked)>1e6) {
                    notice("Choose a finite parameter with |c| <= 1000000.");return;
                }
                state.c=clicked;state.preset=-1;changePlane(1);
                notice("Selected c="+complexText(state.c)+"; fixed eta="+complexText(state.alpha)+". Julia opened; C / B returns to the seeded c portrait.");
            } else {
                picked=clicked;orbitVisible=true;refresh();
                const auto result=diagnostic.qqrhDirectModel->orbit(picked,state.c,state.iterations);
                notice("z="+complexText(picked)+" / "+masterOrbitSummary(result)+
                    " / steps="+std::to_string(result.iterations)+" / endpoint="+complexText(result.endpoint)+". Save orbit exports the finite trace.");
            }
            return;
        }
        if(masterMode()) {
            if(!diagnostic.masterModel || !diagnostic.masterModel->validModel) {
                notice("Model unavailable: "+(diagnostic.masterModel ? diagnostic.masterModel->statusText : std::string("missing model")));
                return;
            }
            if(!frame || frame->generation!=generation.load() || frame->state.family!=20 || frame->state.generation!=state.generation || frame->state.plane!=state.plane) {
                notice("Wait for the current candidate preview before selecting a parameter or orbit.");return;
            }
            const C clicked=world(state.cameras[state.plane],x,y,width,canvasH);
            if(state.plane==0) {
                if(!std::isfinite(clicked.real()) || !std::isfinite(clicked.imag()) || std::abs(clicked)>1e6) {
                    notice("Choose a finite parameter with |c| <= 1000000.");return;
                }
                state.c=clicked;state.preset=-1;
                changePlane(1);
                notice("Selected c="+complexText(state.c)+". Julia opened; C / B: c portrait; A: lambda plane.");
            } else if(state.plane==2) {
                if(!std::isfinite(clicked.real()) || !std::isfinite(clicked.imag()) || std::abs(clicked)>1e6) {
                    notice("Choose a finite complex lambda with |lambda| <= 1000000.");return;
                }
                state.alpha=clicked;state.preset=-1;
                changePlane(1);
                const auto result=diagnostic.masterModel->orbitAtLambda(C{},state.c,state.alpha,state.iterations);
                notice("Selected lambda="+complexText(state.alpha)+"; c="+complexText(state.c)+" / seed 0: "+masterOrbitSummary(result)+". A: lambda plane.");
            } else {
                picked=clicked;orbitVisible=true;refresh();
                const auto result=diagnostic.masterModel->orbitAtLambda(masterPhysicalPicked(),state.c,state.alpha,state.iterations);
                notice("chart="+complexText(picked)+" / z="+complexText(masterPhysicalPicked())+" / "+masterOrbitSummary(result)+
                    " / steps="+std::to_string(result.iterations)+" / endpoint="+complexText(result.endpoint)+". Save orbit exports the finite trace.");
            }
            return;
        }
        if(tripleMode()) {
            if(!diagnostic.tripleModel || !diagnostic.tripleModel->validModel) {
                notice("Model unavailable: "+(diagnostic.tripleModel ? diagnostic.tripleModel->statusText : std::string("missing model")));
                return;
            }
            if(!frame || frame->generation!=generation.load() || frame->state.family!=19 || frame->state.generation!=state.generation || frame->state.plane!=state.plane) {
                notice("Wait for the current moment preview before selecting a parameter or orbit.");return;
            }
            const C clicked=world(state.cameras[state.plane],x,y,width,canvasH);
            if(state.plane==0) {
                if(!std::isfinite(clicked.real()) || !std::isfinite(clicked.imag()) || std::abs(clicked)>1e6) {
                    notice("Choose a finite moment parameter with |c| <= 1000000.");return;
                }
                state.c=clicked;state.preset=-1;
                changePlane(1);
                notice("Selected c="+complexText(state.c)+". Julia plane opened; C / B returns to the atlas.");
            } else {
                picked=clicked;orbitVisible=true;refresh();
                const auto result=diagnostic.tripleModel->orbit(picked,state.c,state.iterations);
                const std::string endpoint=result.atInfinity ? "infinity" : result.endpointRepresentable ? complexText(result.endpoint) : "finite projective point";
                notice("z="+complexText(picked)+" / "+momenttriple::statusName(result.status)+
                    " / steps="+std::to_string(result.iterations)+" / endpoint="+endpoint+". Save orbit exports the finite trace.");
            }
            return;
        }
        if(quotientMode()) {
            if(!diagnostic.quotientModel || !diagnostic.quotientModel->validModel) {
                notice("Model unavailable: "+(diagnostic.quotientModel ? diagnostic.quotientModel->statusText : std::string("missing model")));
                return;
            }
            if(!frame || frame->generation!=generation.load() || frame->state.family!=18 || frame->state.plane!=state.plane) {
                notice("Wait for the current quotient preview before selecting a parameter or orbit.");return;
            }
            const C clicked=world(state.cameras[state.plane],x,y,width,canvasH);
            if(state.plane==0) {
                if(!std::isfinite(clicked.real()) || !std::isfinite(clicked.imag()) || std::abs(clicked)>1e6) {
                    notice("Choose a finite quotient parameter with |c| <= 1000000.");return;
                }
                state.c=clicked;state.preset=-1;
                changePlane(1);
                notice("Selected c="+complexText(state.c)+". Julia plane opened; C / B returns to the atlas.");
            } else {
                picked=clicked;orbitVisible=true;refresh();
                const auto result=diagnostic.quotientModel->orbit(picked,state.c,state.iterations);
                const std::string endpoint=result.atInfinity ? "infinity" : result.endpointRepresentable ? complexText(result.endpoint) : "finite projective point";
                notice("z="+complexText(picked)+" / "+momentquotient::statusName(result.status)+
                    " / steps="+std::to_string(result.iterations)+" / endpoint="+endpoint+". Save orbit exports the finite trace.");
            }
            return;
        }
        if(cubicMode()) {
            if(!diagnostic.cubicModel || !diagnostic.cubicModel->validModel) {
                notice("Model unavailable: "+(diagnostic.cubicModel ? diagnostic.cubicModel->statusText : std::string("missing model")));
                return;
            }
            if(!frame || frame->generation!=generation.load() || frame->state.family!=17 || frame->state.plane!=state.plane) {
                notice("Wait for the current cubic preview before selecting a parameter or orbit.");return;
            }
            const C clicked=world(state.cameras[state.plane],x,y,width,canvasH);
            if(state.plane==0) {
                if(!std::isfinite(clicked.real()) || !std::isfinite(clicked.imag()) || std::abs(clicked)>1e6) {
                    notice("Choose a finite cubic parameter with |c| <= 1000000.");return;
                }
                state.c=clicked;state.preset=-1;
                changePlane(1);
                notice("Selected c="+complexText(state.c)+". Julia plane opened; C / B returns to the atlas.");
            } else {
                picked=clicked;orbitVisible=true;refresh();
                const auto result=diagnostic.cubicModel->orbit(picked,state.c,state.iterations);
                notice("z="+complexText(picked)+" / "+momentcubic::statusName(result.status)+
                    " / steps="+std::to_string(result.iterations)+". Save orbit exports this finite trace.");
            }
            return;
        }
        if(riemannMode() && (!diagnostic.riemannModel || !diagnostic.riemannModel->validModel)) {
            notice("Model unavailable: "+(diagnostic.riemannModel ? diagnostic.riemannModel->statusText : std::string("missing model")));
            return;
        }
        if(projectMode() && (!diagnostic.projectModel || !diagnostic.projectModel->validModel)) {
            notice("Model unavailable: "+(diagnostic.projectModel ? diagnostic.projectModel->statusText : std::string("missing model")));
            return;
        }
        if (spinMode()) {
            if (spinAtlas()) {
                if (!frame || frame->generation!=generation.load()) {
                    notice("Wait for the current atlas preview before selecting a parameter.");return;
                }
                const auto b=spinPlotBox();
                if (!b.has(x,y)) return;
                state.spin.twist=state.spin.kMin+std::clamp(double(x-b.x)/b.w,0.,1.)*(state.spin.kMax-state.spin.kMin);
                state.spin.damping=state.spin.tauMax-std::clamp(double(y-b.y)/b.h,0.,1.)*(state.spin.tauMax-state.spin.tauMin);
                state.preset=-1;
                refresh();request();
                notice("Selected k="+num(state.spin.twist,6)+", tau="+num(state.spin.damping,6)+". E explores its sphere; B returns here.");
            } else notice(state.generation==3 ? "Drag rotates trajectories. B returns to the parameter atlas; N changes initial points."
                                               : "Drag rotates trajectories. E/B adjusts sample count; N changes initial points.");
            return;
        }
        if (sphereMode()) {notice("Drag rotates the cached sample; N draws a new reproducible sample. A shows both hemispheres.");return;}
        if (quantumMode()) {notice("Drag rotates the histories; N draws a new reproducible sample. A compares the ensemble mean.");return;}
        if (state.family == 9) {
            int pw = std::max(1, width / 3), local = x % pw;
            picked = world(state.cameras[0], local, y, pw, canvasH);
            auto e = grh_observable::evaluate(picked.real(), picked.imag());
            if (e.valid)
                notice("a=" + num(e.a, 6) + "   r=" + num(e.r, 9) + "   C=" + complexText(e.observable) +
                       "   symmetry error=" + num(e.symmetryError, 3));
            else
                notice("Observable outside validated numerical window at x=" + num(picked.real(), 4) +
                       ", y=" + num(picked.imag(), 4));
            return;
        }
        picked = world(state.cameras[state.plane], x, y, width, canvasH);
        if (state.family == 11 && state.plane == 0) {
            state.seed = picked;
            state.preset = -1;
        }
        if (state.plane == 1)
            state.c = picked;
        if (state.plane == 2)
            state.alpha = picked;
        if (state.plane)
            state.preset = -1;
        orbitVisible = true;
        refresh();
        if (state.plane)
            request();
        else
            InvalidateRect(window, nullptr, FALSE);
        if (state.family == 11 && state.plane == 0) {
            notice("Orbit traced; parameter seed = " + complexText(picked) + ". C / A opens its parameter slice.");
            return;
        }
        notice(state.plane ? "Parameter selected. E / Explore point opens its dynamical plane."
                           : "Orbit traced from z = " + complexText(picked) +
                                 ". Save orbit exports homogeneous coordinates.");
    }
    void key(WPARAM k) {
        bool ctrl = GetKeyState(VK_CONTROL) < 0;
        if (inCatalogue) {
            if (k == VK_UP)
                selectEntry(selectedEntry - 1);
            else if (k == VK_DOWN)
                selectEntry(selectedEntry + 1);
            else if (k == VK_HOME)
                selectEntry(0);
            else if (k == VK_END)
                selectEntry(int(entries.size()) - 1);
            else if (k == VK_PRIOR)
                selectEntry(selectedEntry - visibleEntries());
            else if (k == VK_NEXT)
                selectEntry(selectedEntry + visibleEntries());
            else if (k == VK_RETURN)
                openEntry();
            else if (k == 'O' && ctrl)
                viewDialog(false);
            else if (k == VK_F11)
                toggleFullscreen();
            else if (k == VK_ESCAPE && full)
                toggleFullscreen();
            return;
        }
        if (k == VK_TAB) {
            showCatalogue();
            return;
        }
        if (k >= '1' && k <= '4') {
            if(masterMode() && labControls) {applyLabPreset(int(k-'1'));return;}
            preset(state, int(k - '1'));
            orbitVisible = false;
            refresh();
            request();
            return;
        }
        if(connectedLambdaMode()) {
            if(k=='F')return;
            if(k=='B' || k=='C' || k=='K') {
                if(k=='B')state.orbitMode=(state.orbitMode+1)%3;
                else if(k=='C')state.spectralOrder=state.spectralOrder%3+1;
                else state.method=(state.method+1)%3;
                state.preset=-1;refresh();request();return;
            }
        }
        if(legacyLambdaMode()) {
            if(k=='F') {setValue(4,(state.method+1)%5);return;}
            if(k=='B') {setValue(5,(state.orbitMode+1)%3);return;}
            if(k=='C') {setValue(6,state.spectralOrder%3);return;}
            if(k=='P' || k=='M') {state.palette=(state.palette+1)%9;request();return;}
            if(k=='A' || k=='J' || k=='E')return;
        }
        if (sampleMode()) {
            if (k=='N') {newSampleHistory();return;}
            if (k==VK_LEFT || k==VK_RIGHT || k==VK_UP || k==VK_DOWN) {
                if (spinAtlas()) {
                    const double dx=(k==VK_LEFT ? -.08 : k==VK_RIGHT ? .08 : 0)*(state.spin.kMax-state.spin.kMin);
                    const double dy=(k==VK_DOWN ? -.08 : k==VK_UP ? .08 : 0)*(state.spin.tauMax-state.spin.tauMin);
                    panSpinAtlas(dx,dy);return;
                }
                if(k==VK_LEFT || k==VK_RIGHT) state.quantum.yaw=std::remainder(state.quantum.yaw+(k==VK_LEFT ? -.1:.1),2*pi);
                else state.quantum.pitch=std::clamp(state.quantum.pitch+(k==VK_UP ? .1:-.1),-pi/2,pi/2);
                request();return;
            }
        }
        if(directQqrhMode()) {
            if(k=='A' || k=='F' || k=='M')return;
            if(k=='K') {state.method=(state.method+1)%3;state.preset=-1;refresh();request();return;}
        }
        if(extraMode() && k=='P') {
            state.palette=(state.palette+1)%3;state.preset=-1;request();return;
        }
        if(hilbertField() && k=='I') {
            hilbertInspectorVisible=!hilbertInspectorVisible;
            InvalidateRect(window,nullptr,FALSE);return;
        }
        if(masterMode()) {
            if(k=='V') {toggleLabControls();return;}
            if(k=='T') {cycleFocusBank();return;}
            if(k=='F') {cycleLab(0);return;}
            if(k=='D') {cycleLab(1);return;}
            if(k=='I') {cycleLab(2);return;}
        }
        if(masterMode() || directQqrhMode()) {
            if(k=='K') {state.spectralOrder=state.spectralOrder%3+1;state.preset=-1;refresh();request();return;}
            if(k=='P') {state.palette=(state.palette+1)%3;request(true);return;}
            if(k==VK_OEM_4 || k==VK_OEM_6) {
                state.hue=std::fmod(state.hue+(k==VK_OEM_4 ? -.025 : .025)+1.,1.);
                request(true);return;
            }
            if(k==VK_OEM_COMMA || k==VK_OEM_PERIOD) {
                state.bands=std::clamp(state.bands+(k==VK_OEM_COMMA ? -.05 : .05),.2,3.);
                request(true);notice("Band frequency="+num(state.bands,2)+". Comma / period adjusts the colour bands.");return;
            }
        }
        switch (k) {
        case 'J':
            if(directAtlasMode()) changePlane(1);
            else if (state.family != 9) changePlane(0);
            break;
        case 'C':
            if(directAtlasMode()) changePlane(0);
            else if (state.family != 9) changePlane(1);
            break;
        case 'A':
            if (state.family != 9) changePlane(2);
            break;
        case 'E':
            explore();
            break;
        case 'B':
            goBack();
            break;
        case 'R':
            resetCamera();
            break;
        case 'H':
            help = !help;
            InvalidateRect(window, nullptr, FALSE);
            break;
        case 'S':
            if (ctrl)
                viewDialog(true);
            else
                saveBmp();
            break;
        case 'O':
            if (ctrl)
                viewDialog(false);
            break;
        case VK_F11:
            toggleFullscreen();
            break;
        case VK_ESCAPE:
            if (hqBusy) {
                ++exportGeneration;
                notice("Cancelling HQ...");
            } else if (help) {
                help = false;
                InvalidateRect(window, nullptr, FALSE);
            } else if (full)
                toggleFullscreen();
            else
                showCatalogue();
            break;
        case VK_ADD:
        case VK_OEM_PLUS:
            zoom(width / 2, canvasH / 2, .8);
            break;
        case VK_SUBTRACT:
        case VK_OEM_MINUS:
            zoom(width / 2, canvasH / 2, 1.25);
            break;
        case VK_LEFT:
            state.cameras[state.plane].x -= state.cameras[state.plane].span * .08;
            request();
            break;
        case VK_RIGHT:
            state.cameras[state.plane].x += state.cameras[state.plane].span * .08;
            request();
            break;
        case VK_UP:
            state.cameras[state.plane].y += state.cameras[state.plane].span * .08;
            request();
            break;
        case VK_DOWN:
            state.cameras[state.plane].y -= state.cameras[state.plane].span * .08;
            request();
            break;
        case VK_HOME:
            resetCamera();
            break;
        }
    }
    void captureClient(const fs::path &path) {
        HDC dc = GetDC(window), mem = CreateCompatibleDC(dc);
        BITMAPINFO bi{};
        bi.bmiHeader.biSize = 40;
        bi.bmiHeader.biWidth = width;
        bi.bmiHeader.biHeight = -height;
        bi.bmiHeader.biPlanes = 1;
        bi.bmiHeader.biBitCount = 32;
        void *pixels = nullptr;
        HBITMAP b = CreateDIBSection(dc, &bi, DIB_RGB_COLORS, &pixels, nullptr, 0);
        auto old = SelectObject(mem, b);
        paint(mem);
        GdiFlush();
        std::vector<uint32_t> data(size_t(width) * height);
        std::memcpy(data.data(), pixels, data.size() * 4);
        writeBmp(path, data, width, height);
        SelectObject(mem, old);
        DeleteObject(b);
        DeleteDC(mem);
        ReleaseDC(window, dc);
    }
    LRESULT handle(UINT msg, WPARAM wp, LPARAM lp) {
        switch (msg) {
        case WM_ERASEBKGND:
            return 1;
        case WM_PAINT: {
            PAINTSTRUCT p;
            HDC dc = BeginPaint(window, &p);
            paint(dc);
            EndPaint(window, &p);
            return 0;
        }
        case WM_SIZE:
            if (wp != SIZE_MINIMIZED) {
                width = LOWORD(lp);
                height = HIWORD(lp);
                canvasH = std::max(sc(120), height - sc(422 + panelExtra()));
                if (inCatalogue) {
                    catalogueScroll =
                        std::clamp(catalogueScroll, 0, std::max(0, int(entries.size()) - visibleEntries()));
                    if (selectedEntry >= catalogueScroll + visibleEntries())
                        catalogueScroll = selectedEntry - visibleEntries() + 1;
                }
                request();
            }
            return 0;
        case WM_GETMINMAXINFO: {
            auto p = reinterpret_cast<MINMAXINFO *>(lp);
            p->ptMinTrackSize = {sc(1000), sc(645)};
            return 0;
        }
        case WM_LBUTTONDBLCLK:
            if (inCatalogue)
                catalogueClick(GET_X_LPARAM(lp), GET_Y_LPARAM(lp), true);
            return 0;
        case WM_LBUTTONDOWN:
            click(GET_X_LPARAM(lp), GET_Y_LPARAM(lp));
            return 0;
        case WM_MOUSEMOVE:
            if (active >= 0)
                moveSlider(GET_X_LPARAM(lp));
            else if (dragging) {
                if (spinAtlas()) {
                    state.spin.kMin=dragSpin.kMin;state.spin.kMax=dragSpin.kMax;
                    state.spin.tauMin=dragSpin.tauMin;state.spin.tauMax=dragSpin.tauMax;
                    panSpinAtlas(-(GET_X_LPARAM(lp)-dragX)*(dragSpin.kMax-dragSpin.kMin)/dragSpinPlot.w,
                                  (GET_Y_LPARAM(lp)-dragY)*(dragSpin.tauMax-dragSpin.tauMin)/dragSpinPlot.h);
                    return 0;
                }
                if(legacyLambdaMode()) {
                    const auto b=legacyLambdaBox();
                    auto &cam=state.cameras[2];
                    cam.x=dragCamera.x-(GET_X_LPARAM(lp)-dragX)*dragCamera.span/b.w;
                    cam.y=dragCamera.y+(GET_Y_LPARAM(lp)-dragY)*dragCamera.span/b.h;
                    request();return 0;
                }
                if (sampleMode()) {
                    state.quantum.yaw=std::remainder(dragCamera.x+(GET_X_LPARAM(lp)-dragX)*pi/std::max(100,width),2*pi);
                    state.quantum.pitch=std::clamp(dragCamera.y+(GET_Y_LPARAM(lp)-dragY)*pi/std::max(100,canvasH),-pi/2,pi/2);
                    request();return 0;
                }
                auto &c = state.cameras[state.plane];
                int mapW = state.family == 9 ? std::max(1, width / 3) : width;
                c.x = dragCamera.x - (GET_X_LPARAM(lp) - dragX) * dragCamera.span / mapW;
                c.y = dragCamera.y + (GET_Y_LPARAM(lp) - dragY) * dragCamera.span / mapW;
                if (state.family == 9) state.cameras[1] = state.cameras[2] = c;
                request();
            }
            return 0;
        case WM_LBUTTONUP:
            active = -1;
            dragging = false;
            ReleaseCapture();
            finishInteraction();
            InvalidateRect(window, nullptr, FALSE);
            return 0;
        case WM_CAPTURECHANGED: {
            const bool wasInteracting=active>=0 || dragging;
            active = -1;
            dragging = false;
            if(wasInteracting)finishInteraction();
            return 0;
        }
        case WM_RBUTTONDOWN:
            rightClick(GET_X_LPARAM(lp), GET_Y_LPARAM(lp));
            return 0;
        case WM_MOUSEWHEEL: {
            POINT p{GET_X_LPARAM(lp), GET_Y_LPARAM(lp)};
            ScreenToClient(window, &p);
            wheelInput=true;
            try {wheel(p.x, p.y, GET_WHEEL_DELTA_WPARAM(wp));}
            catch(...) {wheelInput=false;throw;}
            wheelInput=false;
            return 0;
        }
        case WM_TIMER:
            if(wp==INTERACTION_TIMER) {interactionTick();return 0;}
            break;
        case WM_KEYDOWN:
            key(wp);
            return 0;
        case FRAME: {
            std::unique_ptr<Frame> f(reinterpret_cast<Frame *>(lp));
            if (f->generation == generation.load()) {
                backend = f->backend;
                frame = std::move(f);
                if(frame->stride==1 && frame->state.family==20) {
                    fullOrbitState=frame->state;
                    fullOrbitWidth=frame->targetW;fullOrbitHeight=frame->targetH;
                }
                InvalidateRect(window, nullptr, FALSE);
                if (!options.uiSmoke.empty() && frame->stride == 1) {
                    UpdateWindow(window);
                    captureClient(options.uiSmoke);
                    std::cout << "UI_SMOKE_PASS=" << utf8(options.uiSmoke.wstring()) << std::endl;
                    PostMessageW(window, WM_CLOSE, 0, 0);
                }
            }
            return 0;
        }
        case RENDER_DONE:
            if(previewInFlight && previewGeneration==static_cast<uint64_t>(wp)) {
                previewInFlight=false;
                interactionTick();
            }
            return 0;
        case GPU_READY:
            notice(wp?"GPU shaders ready. The selected backend can now use them.":
                "GPU preparation failed. CPU rendering remains available.");
            if(wp&&!options.cpu&&!rhPreferCpu(options,state)&&!inCatalogue)request();
            return 0;
        case NOTICE: {
            std::unique_ptr<std::string> s(reinterpret_cast<std::string *>(lp));
            if (wp == 1)
                hqBusy = false;
            notice(*s);
            return 0;
        }
        case WM_CLOSE:
            DestroyWindow(window);
            return 0;
        case WM_DESTROY: {
            KillTimer(window,INTERACTION_TIMER);
            releasePaintBuffer();
            std::lock_guard<std::mutex> lock(mutex);
            closing = true;
            ++generation;
            ++exportGeneration;
            cv.notify_one();
            PostQuitMessage(0);
            return 0;
        }
        }
        return DefWindowProcW(window, msg, wp, lp);
    }
    static LRESULT CALLBACK proc(HWND h, UINT m, WPARAM w, LPARAM l) {
        App *app = reinterpret_cast<App *>(GetWindowLongPtrW(h, GWLP_USERDATA));
        if (m == WM_NCCREATE) {
            app = static_cast<App *>(reinterpret_cast<CREATESTRUCTW *>(l)->lpCreateParams);
            app->window = h;
            SetWindowLongPtrW(h, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(app));
        }
        if (app) {
            try {
                return app->handle(m, w, l);
            } catch (const std::exception &e) {
                app->notice("Error: " + std::string(e.what()));
                return 0;
            }
        }
        return DefWindowProcW(h, m, w, l);
    }

  public:
    App(Options o, fs::path p) : options(std::move(o)), root(std::move(p)) {
        std::vector<std::string> warnings;
        entries = catalogue::load(root / L"renderers", &warnings);
        for (auto &warning : warnings)
            std::cerr << "CATALOGUE_WARNING=" << warning << std::endl;
        for (auto &e : entries) {
            auto initial=defaults(e.family,e.generation);
            if(e.family==20)focus_presets::apply(initial,1,0);
            entryStates.push_back(initial);
        }
        state = entryStates[0];
        if(options.renderer.empty()&&options.view.empty()&&entries.size()>std::size_t(rhfresh::OriginalCandidateCount)) {
            selectedEntry=rhfresh::CandidateCount;
            catalogueScroll=std::min(selectedEntry,std::max(0,int(entries.size())-visibleEntries()));
            state=entryStates[selectedEntry];
        }
        if (!options.renderer.empty()) {
            auto found = std::find_if(entries.begin(), entries.end(),
                                      [&](const catalogue::Entry &e) { return e.id == options.renderer; });
            if (found == entries.end())
                throw std::runtime_error("Unknown renderer ID: " + options.renderer);
            selectedEntry = activeEntry = int(found - entries.begin());
            state = entryStates[selectedEntry];
            inCatalogue = false;
        }
        if (!options.view.empty()) {
            state = readView(options.view);
            identifyState();
            inCatalogue = false;
        }
        refresh();
    }
    int run() {
        SetProcessDPIAware();
        HDC screen = GetDC(nullptr);
        dpiScale = GetDeviceCaps(screen, LOGPIXELSX) / 96.f;
        ReleaseDC(nullptr, screen);
        font =
            CreateFontW(-sc(13), 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS,
                        CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, FIXED_PITCH | FF_MODERN, L"Consolas");
        smallFont =
            CreateFontW(-sc(12), 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS,
                        CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY, FIXED_PITCH | FF_MODERN, L"Consolas");
        titleFont = CreateFontW(-sc(25), 0, 0, 0, FW_SEMIBOLD, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
                                OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
                                FIXED_PITCH | FF_MODERN, L"Consolas");
        WNDCLASSW wc{};
        wc.style = CS_DBLCLKS;
        wc.lpfnWndProc = proc;
        wc.hInstance = GetModuleHandleW(nullptr);
        wc.hCursor = LoadCursorW(nullptr, IDC_ARROW);
        wc.hIcon = LoadIconW(nullptr, IDI_APPLICATION);
        wc.lpszClassName = L"AmpNativeWindow";
        RegisterClassW(&wc);
        RECT work;
        SystemParametersInfoW(SPI_GETWORKAREA, 0, &work, 0);
        int w = std::min(sc(1280), int(work.right - work.left) - sc(28)),
            h = std::min(sc(976), int(work.bottom - work.top) - sc(38));
        window = CreateWindowExW(0, wc.lpszClassName, title, WS_OVERLAPPEDWINDOW, CW_USEDEFAULT,
                                 CW_USEDEFAULT, w, h, nullptr, nullptr, wc.hInstance, this);
        if (!window)
            throw std::runtime_error("Cannot create native window");
        HMODULE dwm = LoadLibraryW(L"dwmapi.dll");
        if (dwm) {
            using DwmAttr = HRESULT(WINAPI *)(HWND, DWORD, LPCVOID, DWORD);
            auto fn = reinterpret_cast<DwmAttr>(GetProcAddress(dwm, "DwmSetWindowAttribute"));
            BOOL dark = TRUE;
            if (fn)
                fn(window, 20, &dark, sizeof(dark));
            FreeLibrary(dwm);
        }
        ShowWindow(window, SW_SHOWNORMAL);
        UpdateWindow(window);
        worker = std::thread([this] { workerLoop(); });
        request();
        MSG msg;
        while (GetMessageW(&msg, nullptr, 0, 0) > 0) {
            TranslateMessage(&msg);
            DispatchMessageW(&msg);
        }
        int exitCode = int(msg.wParam);
        if (worker.joinable())
            worker.join();
        while (PeekMessageW(&msg, nullptr, FRAME, RENDER_DONE, PM_REMOVE)) {
            if (msg.message == FRAME)
                delete reinterpret_cast<Frame *>(msg.lParam);
            if (msg.message == NOTICE)
                delete reinterpret_cast<std::string *>(msg.lParam);
        }
        releasePaintBuffer();
        DeleteObject(font);
        DeleteObject(smallFont);
        DeleteObject(titleFont);
        return exitCode;
    }
};

static int compareBackends() {
    Renderer gpu, cpu;
    cpu.init(true);
    if (!gpu.init(false))
        throw std::runtime_error("GPU comparison unavailable: " + gpu.diagnostic);
    std::cout << "COMPARE_BACKEND=" << gpu.backend << std::endl;
    std::vector<State> cases;
    for (int plane = 0; plane < 3; plane++)
        for (int palette = 0; palette < 3; palette++) {
            State s;
            s.plane = plane;
            s.palette = palette;
            s.iterations = 96;
            cases.push_back(s);
        }
    for (int family = 1; family <= 5; family++)
        for (int gen = 1; gen <= 4; gen++)
            for (int method = 0; method < 2; method++) {
                auto s = defaults(family, gen);
                s.method = method;
                s.iterations = 96;
                cases.push_back(s);
                if (gen == 2)
                    for (int plane = 1; plane < 3; plane++) {
                        auto p = s;
                        p.plane = plane;
                        cases.push_back(p);
                    }
            }
    for (int method = 0; method < 2; method++) {
        for (int palette = 1; palette <= 2; palette++) {
            auto s = defaults(3, 3);
            s.method = method;
            s.palette = palette;
            s.iterations = 96;
            cases.push_back(s);
        }
        auto s = defaults(3, 4);
        s.c = {0, 0};
        s.method = method;
        s.iterations = 96;
        cases.push_back(s);
    }
    for (int presetIndex = 1; presetIndex <= 2; presetIndex++) {
        auto s = defaults(3, 2);
        preset(s, presetIndex);
        s.iterations = 96;
        cases.push_back(s);
    }
    for (int method = 0; method < 2; method++) {
        for (int palette = 1; palette <= 2; palette++) {
            auto s = defaults(4, 3);
            s.method = method;
            s.palette = palette;
            s.iterations = 96;
            cases.push_back(s);
        }
        auto s = defaults(4, 1);
        s.c = {0, 0};
        s.method = method;
        s.iterations = 256;
        cases.push_back(s);
    }
    for (int p = 1; p <= 2; p++) {
        auto s = defaults(4, 2);
        preset(s, p);
        s.iterations = 96;
        cases.push_back(s);
    }
    int passed = 0;
    for (const auto &s : cases) {
        auto r = scene(s);
        int w = s.family ? 96 : 160, h = s.family ? 64 : 96;
        std::vector<uint32_t> a, b;
        gpu.render(r, w, h, double(h) / w, 1, a);
        if (!gpu.gpu)
            throw std::runtime_error("Unexpected GPU fallback: " + gpu.diagnostic);
        cpu.render(r, w, h, double(h) / w, 1, b);
        double sum = 0;
        int severe = 0;
        for (size_t i = 0; i < a.size(); i++) {
            int largest = 0;
            for (int k = 0; k < 3; k++) {
                int d = std::abs(int((a[i] >> (8 * k)) & 255) - int((b[i] >> (8 * k)) & 255));
                sum += d;
                largest = std::max(largest, d);
            }
            if (largest > 30)
                severe++;
        }
        double mean = sum / (3 * a.size()), fraction = severe / double(a.size());
        std::cout << "family=" << s.family << " generation=" << s.generation << " method=" << s.method
                  << " plane=" << s.plane << " palette=" << s.palette << " mean_channel_error=" << mean
                  << " large_difference_fraction=" << fraction << std::endl;
        if (mean > 12 || fraction > .15)
            throw std::runtime_error("CPU/GPU comparison outside tolerance");
        passed++;
    }
    std::cout << "PASS: " << passed << " CPU / compute shader image comparisons\n";
    return 0;
}

// Batch rendering runs only in the user's native executable.
// The DIB owns a top-down contact sheet so GDI labels and renderer pixels share coordinates.
struct FiveSheet {
    HDC dc=nullptr;
    HBITMAP bitmap=nullptr;
    HFONT heading=nullptr,body=nullptr;
    HGDIOBJ previousBitmap=nullptr,previousFont=nullptr;
    uint32_t *pixels=nullptr;
    int width=0,height=0;
    FiveSheet(int w,int h,int fontSize):width(w),height(h) {
        dc=CreateCompatibleDC(nullptr);
        BITMAPINFO bi{};
        bi.bmiHeader.biSize=sizeof(BITMAPINFOHEADER);
        bi.bmiHeader.biWidth=w;bi.bmiHeader.biHeight=-h;
        bi.bmiHeader.biPlanes=1;bi.bmiHeader.biBitCount=32;
        bi.bmiHeader.biCompression=BI_RGB;
        void *bits=nullptr;
        if(dc) bitmap=CreateDIBSection(dc,&bi,DIB_RGB_COLORS,&bits,nullptr,0);
        if(!dc || !bitmap || !bits) {release();throw std::runtime_error("Cannot allocate five-map contact sheet");}
        pixels=static_cast<uint32_t*>(bits);
        previousBitmap=SelectObject(dc,bitmap);
        heading=CreateFontW(-fontSize,0,0,0,FW_SEMIBOLD,FALSE,FALSE,FALSE,DEFAULT_CHARSET,
            OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,ANTIALIASED_QUALITY,DEFAULT_PITCH,L"Segoe UI");
        body=CreateFontW(-std::max(12,fontSize*3/4),0,0,0,FW_NORMAL,FALSE,FALSE,FALSE,DEFAULT_CHARSET,
            OUT_DEFAULT_PRECIS,CLIP_DEFAULT_PRECIS,ANTIALIASED_QUALITY,DEFAULT_PITCH,L"Segoe UI");
        if(!heading || !body) {release();throw std::runtime_error("Cannot create contact-sheet fonts");}
        previousFont=SelectObject(dc,heading);
        std::fill(pixels,pixels+size_t(w)*h,0x080b14u);
        SetBkMode(dc,TRANSPARENT);
    }
    FiveSheet(const FiveSheet&)=delete;
    FiveSheet &operator=(const FiveSheet&)=delete;
    ~FiveSheet(){release();}
    void release() {
        if(dc && previousFont) SelectObject(dc,previousFont);
        if(dc && previousBitmap) SelectObject(dc,previousBitmap);
        if(heading) DeleteObject(heading);
        if(body) DeleteObject(body);
        if(bitmap) DeleteObject(bitmap);
        if(dc) DeleteDC(dc);
        dc=nullptr;bitmap=nullptr;heading=nullptr;body=nullptr;pixels=nullptr;
    }
    void text(int x,int y,int w,int h,const std::string &s,bool large=false) {
        SelectObject(dc,large?heading:body);
        SetTextColor(dc,large?RGB(239,229,255):RGB(190,205,229));
        RECT box{x,y,x+w,y+h};const auto value=wide(s);
        DrawTextW(dc,value.c_str(),int(value.size()),&box,DT_LEFT|DT_WORDBREAK|DT_NOPREFIX);
        GdiFlush();
    }
    void copy(int x,int y,int w,int h,const std::vector<uint32_t>&data) {
        GdiFlush();
        for(int row=0;row<h;++row)
            std::copy_n(data.data()+size_t(row)*w,w,pixels+size_t(y+row)*width+x);
    }
    std::vector<uint32_t> data() {GdiFlush();return {pixels,pixels+size_t(width)*height};}
};
static int plotFive(const Options &o) {
    // Keep the overview readable even for small individual images, and bound allocations.
    const int tileW=std::max(o.width,480),imageH=o.height;
    const int gap=16,titleH=58,footerH=38,topH=82;
    const int sheetW=3*tileW+4*gap;
    const int tileH=std::max(titleH+imageH+footerH,512);
    const int sheetH=topH+2*tileH+3*gap;
    if(uint64_t(sheetW)*sheetH>80000000ULL)
        throw std::runtime_error("--plot-five contact sheet exceeds 80 million pixels; reduce --width or --height");
    fs::create_directories(o.plotFive);
    const std::wstring prefix=L"five_projects_"+wide(stamp());
    fs::path folder=o.plotFive/prefix;
    for(unsigned attempt=1;fs::exists(folder);++attempt)
        folder=o.plotFive/(prefix+L"_"+std::to_wstring(attempt));
    if(!fs::create_directory(folder)) throw std::runtime_error("Cannot create batch output folder");
    std::ofstream index(folder/L"INDEX.txt");
    if(!index) throw std::runtime_error("Cannot create five-map index");
    index << "FRAKTALFAKTUR 0.24.0 / five project maps\nIndividual images: " << o.width << 'x' << o.height
          << "\nOne sampled complex start per pixel; palette 0; initial preset of each project.\n"
          << "The GRH picture uses a finite polynomial approximation, not the exact completed-L dynamics.\n"
          << "Dark or cutoff pixels remain unresolved; no escape or conjecture proof follows from these plots.\n\n";
    FiveSheet sheet(sheetW,sheetH,std::clamp(tileW/36,18,32));
    sheet.text(gap,13,sheetW-2*gap,42,"FRAKTALFAKTUR / five project iterations",true);
    sheet.text(gap,53,sheetW-2*gap,28,"GRH    PI    NSE    COLLATZ    P_c(d,q)   |   complex initial-state planes / finite numerical diagnostics");
    static const char *tags[]={"GRH","PI","NSE","COLLATZ","PCD"};
    Renderer engine;engine.init(o.cpu,o.threads);engine.setGpuTuning(o.gpuTuning);
    for(int i=0;i<5;++i) {
        State s=defaults(15,i+1);
        auto r=scene(s);
        if(!r.projectModel || !r.projectModel->validModel)
            throw std::runtime_error(std::string(tags[i])+" model is unavailable: "+(r.projectModel?r.projectModel->statusText:"missing model"));
        const auto filename=wide("0"+std::to_string(i+1)+"_"+tags[i]+".bmp");
        const auto path=folder/filename;
        const auto begin=Clock::now();
        std::cout << "PLOT=" << i+1 << "/5 PROJECT=" << tags[i] << " SIZE=" << o.width << 'x' << o.height
                  << " ITERATIONS=" << s.iterations << " BACKEND=" << engine.currentBackend(r,o.width) << std::endl;
        std::vector<uint32_t> data;
        if(!engine.render(r,o.width,o.height,double(o.height)/o.width,1,data))
            throw std::runtime_error(std::string(tags[i])+" batch render did not complete");
        writeBmp(path,data,o.width,o.height);
        metadata(path,s,o.width,o.height,1,engine.currentBackend(r,o.width),nullptr,nullptr,nullptr,r.projectModel.get());
        const int x=gap+(i%3)*(tileW+gap),y=topH+gap+(i/3)*(tileH+gap);
        sheet.text(x,y,tileW,titleH-4,projectTitle(i+1),true);
        sheet.copy(x+(tileW-o.width)/2,y+titleH,o.width,o.height,data);
        const auto &cam=s.cameras[0];
        const auto caption="centre "+complexText(C(cam.x,cam.y))+" | span "+num(cam.span,3)+" | steps "+std::to_string(s.iterations);
        sheet.text(x,y+titleH+imageH+7,tileW,footerH-7,caption);
        index << utf8(filename) << "\n" << r.projectModel->description << "\n" << caption << "\n\n";
        std::cout << "SAVED=" << utf8(fs::absolute(path).wstring()) << "\nPLOT_MS="
                  << std::chrono::duration<double,std::milli>(Clock::now()-begin).count() << std::endl;
    }
    const int lx=gap+2*(tileW+gap),ly=topH+gap+tileH+gap;
    sheet.text(lx,ly,tileW,titleH,"How to read these five fields",true);
    sheet.text(lx,ly+titleH,tileW,tileH-titleH,
        "GRH: rational dynamics of a finite Taylor model of the completed character pair.\n\n"
        "PI: Newton basins from actual binary pi pairs with the analytic drift removed.\n\n"
        "NSE: Newton basins of a diffusing finite spectral-energy polynomial.\n\n"
        "COLLATZ: two steps of the entire shortcut extension.\n\n"
        "P_c(d,q): Newton basins of the q=2 GP-4 critical-dimension quintic.\n\n"
        "Colours show detected roots or cycles and finite orbit structure. Numerical cutoffs remain unclassified.\n\n"
        "Each BMP has a .txt record and a reloadable .view file. Open the corresponding catalogue entry to pan, zoom and trace orbits.");
    const auto contact=folder/L"index.bmp";
    writeBmp(contact,sheet.data(),sheetW,sheetH);
    index << "Contact sheet: index.bmp\nNative user-side rendering completed.\n";
    if(!index) throw std::runtime_error("Five-map index write failed");
    index.close();
    std::cout << "PLOT_FIVE=PASS\nCONTACT_SHEET=" << utf8(fs::absolute(contact).wstring())
              << "\nOUTPUT_DIRECTORY=" << utf8(fs::absolute(folder).wstring()) << std::endl;
    const auto opened=ShellExecuteW(nullptr,L"open",fs::absolute(contact).c_str(),nullptr,nullptr,SW_SHOWNORMAL);
    if(reinterpret_cast<INT_PTR>(opened)<=32)
        std::cout << "OPEN_CONTACT_SHEET=MANUAL (open index.bmp in the output folder)" << std::endl;
    return 0;
}

int wmain(int argc, wchar_t **argv) {
    try {
        Options o;
        for (int i = 1; i < argc; i++) {
            std::wstring a = argv[i];
            auto arg = [&]() -> std::wstring {
                if (i + 1 >= argc)
                    throw std::runtime_error("Missing argument");
                return argv[++i];
            };
            if (a == L"--cpu")
                o.cpu = true;
            else if (a == L"--self-test")
                o.test = true;
            else if (a == L"--gpu-self-test")
                o.gpuTest = true;
            else if (a == L"--benchmark")
                o.benchmark = true;
            else if (a == L"--engine-probe")
                o.engineProbe = true;
            else if (a == L"--probe-full") {
                o.engineProbe=true; o.probeOptions.full=true;
            } else if (a == L"--probe-repeats")
                o.probeOptions.repeats=std::stoi(arg());
            else if (a == L"--probe-seconds")
                o.probeOptions.maxSeconds=std::stod(arg());
            else if (a == L"--probe-iterations")
                o.probeOptions.iterations=std::stoi(arg());
            else if (a == L"--gpu-schedule") {
                const auto mode=arg();
                if(mode==L"legacy")o.gpuTuning.schedule=rhgpu::Schedule::Legacy;
                else if(mode==L"adaptive")o.gpuTuning.schedule=rhgpu::Schedule::Adaptive;
                else throw std::runtime_error("--gpu-schedule must be legacy or adaptive");
            } else if (a == L"--gpu-batch-ms")
                o.gpuTuning.targetBatchMs=std::stod(arg());
            else if (a == L"--gpu-wait") {
                const auto mode=arg();
                if(mode==L"responsive")o.gpuTuning.waitPolicy=rhgpu::WaitPolicy::Responsive;
                else if(mode==L"legacy-sleep")o.gpuTuning.waitPolicy=rhgpu::WaitPolicy::LegacySleep;
                else throw std::runtime_error("--gpu-wait must be responsive or legacy-sleep");
            } else if (a == L"--gpu-wait-ms")
                o.gpuTuning.waitSliceMs=std::stod(arg());
            else if (a == L"--gpu-min-pixels")
                o.gpuTuning.minimumBatchPixels=unsigned(std::stoul(arg()));
            else if (a == L"--gpu-batch-pixels")
                o.gpuTuning.batchPixels=std::stoull(arg());
            else if (a == L"--stable-curvature")
                o.gpuTuning.stableCurvature=true;
            else if (a == L"--rh-backend") {
                o.rhBackend=utf8(arg());
                if(o.rhBackend!="auto"&&o.rhBackend!="cpu"&&o.rhBackend!="gpu")
                    throw std::runtime_error("--rh-backend must be auto, cpu or gpu");
            }
            else if (a == L"--no-colour-cache")
                o.gpuTuning.metricCache=false;
            else if (a == L"--interaction") {
                o.interaction=utf8(arg());
                if(o.interaction!="coarse"&&o.interaction!="coalesced"&&o.interaction!="latest")
                    throw std::runtime_error("--interaction must be coarse, coalesced or latest");
            } else if (a == L"--preview-device") {
                o.previewDevice=utf8(arg());
                if(o.previewDevice!="auto"&&o.previewDevice!="cpu"&&o.previewDevice!="gpu")
                    throw std::runtime_error("--preview-device must be auto, cpu or gpu");
            }
            else if (a == L"--compare-backends")
                throw std::runtime_error("Use the registered RH candidate renderers.");
            else if (a == L"--list")
                o.list = true;
            else if (a == L"--renderer")
                o.renderer = utf8(arg());
            else if (a == L"--threads") {
                int value = std::stoi(arg());
                if (value < 1 || value > 64)
                    throw std::runtime_error("--threads must be 1..64");
                o.threads = unsigned(value);
            } else if (a == L"--refine") {
                o.refine = utf8(arg());
                if (o.refine != "default" && o.refine != "fast" && o.refine != "direct")
                    throw std::runtime_error("--refine must be default, fast, or direct");
            } else if (a == L"--view")
                o.view = arg();
            else if (a == L"--plot-five")
                throw std::runtime_error("Use --renderer with a registered RH candidate ID.");
            else if (a == L"--render")
                o.render = arg();
            else if (a == L"--ui-smoke")
                o.uiSmoke = arg();
            else if (a == L"--width") {
                o.width = std::stoi(arg());o.probeOptions.width=o.width;
            } else if (a == L"--height") {
                o.height = std::stoi(arg());o.probeOptions.height=o.height;
            }
            else if (a == L"--help") {
                std::cout << "FRAKTALFAKTUR RH Lab GPU 1.8.6 / run.cmd [options]\n--list  "
                             "--renderer curvature  --cpu  --threads 32  "
                             "--refine default|fast|direct  --view file.view  --self-test  --gpu-self-test  --benchmark "
                             "\n--render file.bmp --width 1920 "
                             "--height 1080\n--ui-smoke "
                             "screenshot.bmp\n--engine-probe [--probe-full] [--probe-repeats 3] [--probe-seconds 90]\n"
                             "--gpu-schedule legacy|adaptive --gpu-batch-ms 4 --gpu-min-pixels 2048\n"
                             "--gpu-wait responsive|legacy-sleep --gpu-wait-ms 0.35\n"
                             "--preview-device auto|cpu|gpu\n"
                             "--rh-backend auto|cpu|gpu --stable-curvature --gpu-batch-pixels 32768\n"
                             "--interaction coarse|coalesced|latest --no-colour-cache\n";
                return 0;
            } else
                throw std::runtime_error("Unknown argument: " + utf8(a));
        }
        if (o.test) {
            const bool mathPass=rhfresh::runSelfTests(std::cout);
            const bool classicPass=fresh_classic_tests::run(std::cout);
            const bool labPass=lab_tests::run(std::cout);
            const bool focusPass=focus_tests::run(std::cout);
            const bool cpuOptimisationPass=cpu_optimisation_tests::run(std::cout);
            const bool workersPass=render_worker_tests::run(std::cout);
            const bool controlsPass=engine_controls_tests::run(std::cout);
            const bool curvaturePass=stable_curvature_tests::run(std::cout);
            const bool branchPass=rh_branch_extension_tests::run(std::cout);
            const bool newStructurePass=new_structure_tests::run(std::cout);
            const bool extraPass=frk_extra_tests::run(std::cout);
            const bool hilbertPass=qqrh_hilbert_tests::run(std::cout);
            return mathPass && classicPass && labPass && focusPass && cpuOptimisationPass && workersPass && controlsPass && curvaturePass && branchPass && newStructurePass && extraPass && hilbertPass ? 0 : 1;
        }
        if (o.gpuTest) {
            const int formulas=gpu_tests::run(std::cout);
            if(formulas!=0)return formulas;
            const int performance=performance_gpu_tests::run(std::cout);
            if(performance!=0)return performance;
            const int curvature=stable_curvature_gpu_tests::run(std::cout);
            if(curvature!=0)return curvature;
            const int branch=rh_branch_gpu_tests::run(std::cout);
            if(branch!=0)return branch;
            const int newStructureGpu=new_structure_gpu_tests::run(std::cout);
            if(newStructureGpu!=0)return newStructureGpu;
            Com<ID3DBlob> code,error;
            HRESULT compiled=D3DCompile(frkExtraShader,std::strlen(frkExtraShader),"frk_extra_shader",
                                       nullptr,nullptr,"main","cs_5_0",D3DCOMPILE_OPTIMIZATION_LEVEL3,
                                       0,code.put(),error.put());
            if(FAILED(compiled)) {
                std::cout<<"EXTRA_GPU_SHADER=FAIL "<<(error.p ? std::string(
                    static_cast<char*>(error->GetBufferPointer()),error->GetBufferSize()) : "unknown error")<<'\n';
                return 1;
            }
            std::cout<<"EXTRA_GPU_SHADER=PASS BYTECODE_BYTES="<<code->GetBufferSize()<<'\n';
            return qqrh_hilbert_gpu_tests::run(std::cout);
        }
        if (o.compare)
            return compareBackends();
        if(!rhgpu::validTuning(o.gpuTuning))throw std::runtime_error("Invalid GPU tuning");
        if (o.width < 32 || o.height < 32 || o.width > 8192 || o.height > 8192)
            throw std::runtime_error("Render dimensions must be 32..8192");
        wchar_t exe[32768];
        GetModuleFileNameW(nullptr, exe, 32768);
        fs::path root = fs::path(exe).parent_path();
        if(o.engineProbe) {
            if(!o.view.empty())o.probeOptions.customState=readView(o.view);
            return engine_probe::run(std::cout,root,o.probeOptions);
        }
        if(o.benchmark) {
            const bool forceCpu=o.cpu||o.rhBackend!="gpu";
            return performance_benchmark::run(std::cout,forceCpu,forceCpu?rhCpuThreads(o):o.threads,
                o.width,o.height,root/L"performance.csv",o.gpuTuning);
        }
        auto available = catalogue::load(root / L"renderers");
        if (o.list) {
            for (auto &e : available)
                std::cout << e.id << " | " << e.title << " | " << e.description << "\n";
            return 0;
        }
        if(!o.plotFive.empty()) {
            if(!o.render.empty() || !o.renderer.empty() || !o.view.empty() || !o.uiSmoke.empty())
                throw std::runtime_error("--plot-five cannot be combined with --render, --renderer, --view or --ui-smoke");
            return plotFive(o);
        }
        if (!o.render.empty()) {
            State s=defaults(20,1);
            if (!o.renderer.empty()) {
                auto e = std::find_if(available.begin(), available.end(),
                                      [&](const catalogue::Entry &e) { return e.id == o.renderer; });
                if (e == available.end())
                    throw std::runtime_error("Unknown renderer ID: " + o.renderer);
                s = defaults(e->family, e->generation);
            }
            if(s.family==20)focus_presets::apply(s,1,0);
            if (!o.view.empty())
                s = readView(o.view);
            auto r = scene(s);
            if(s.family==24 && (!r.qqrhDirectModel || !r.qqrhDirectModel->validModel))
                throw std::runtime_error("Direct QQRH model unavailable: "+(r.qqrhDirectModel ? r.qqrhDirectModel->statusText : std::string("missing model")));
            if(s.family==20 && (!r.masterModel || !r.masterModel->validModel))
                throw std::runtime_error("RH fresh candidate model unavailable: "+(r.masterModel ? r.masterModel->statusText : std::string("missing model")));
            if(s.family==19 && (!r.tripleModel || !r.tripleModel->validModel))
                throw std::runtime_error("RH moment model unavailable: "+(r.tripleModel ? r.tripleModel->statusText : std::string("missing model")));
            if(s.family==18 && (!r.quotientModel || !r.quotientModel->validModel))
                throw std::runtime_error("RH quotient model unavailable: "+(r.quotientModel ? r.quotientModel->statusText : std::string("missing model")));
            if(s.family==17 && (!r.cubicModel || !r.cubicModel->validModel))
                throw std::runtime_error("RH cubic model unavailable: "+(r.cubicModel ? r.cubicModel->statusText : std::string("missing model")));
            if(s.family==16 && (!r.riemannModel || !r.riemannModel->validModel))
                throw std::runtime_error("RH model unavailable: "+(r.riemannModel ? r.riemannModel->statusText : std::string("missing model")));
            if(s.family==15 && (!r.projectModel || !r.projectModel->validModel))
                throw std::runtime_error("Project model unavailable: "+(r.projectModel ? r.projectModel->statusText : std::string("missing model")));
            Renderer engine;
            const bool forceCpu=o.cpu||rhPreferCpu(o,s);
            engine.init(forceCpu,forceCpu&&s.family==20?rhCpuThreads(o):o.threads);
            engine.setGpuTuning(o.gpuTuning);
            std::vector<uint32_t> data;
            auto start = Clock::now();
            if(!engine.render(r, o.width, o.height, double(o.height) / o.width, 1, data))
                throw std::runtime_error("Render did not complete");
            std::cout << "BACKEND=" << engine.currentBackend(r, o.width)
                      << "\nDETAIL=" << (s.family==20 ? engine.rhGpuDiagnostic() : engine.diagnostic) << std::endl;
            writeBmp(o.render, data, o.width, o.height);
            const auto qs=engine.quantumStats();
            const auto ss=engine.sphericalStats();
            const auto ts=engine.spinStats();
            metadata(o.render, s, o.width, o.height, 1, engine.currentBackend(r, o.width), &qs, &ss, &ts, r.projectModel.get(), r.riemannModel.get(), r.cubicModel.get(), r.quotientModel.get(), r.tripleModel.get(), r.masterModel.get(), r.qqrhDirectModel.get());
            std::cout << "RENDER_MS="
                      << std::chrono::duration<double, std::milli>(Clock::now() - start).count()
                      << "\nSAVED=" << utf8(o.render.wstring()) << std::endl;
            return 0;
        }
        App app(o, root);
        return app.run();
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}
