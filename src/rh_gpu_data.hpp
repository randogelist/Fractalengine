#pragma once

#include "amp_core.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

// Host layout shared with the RH compute shader. One row of the immutable
// RGBA32F texture is one polynomial in x. Row0 also holds the opt-in
// curvature P/H/J polynomials in y/z/w; default evaluation still reads x.
// The coefficient values and radial derivatives come from the same finite
// rounded Taylor table used by the CPU, before rounding to float32 here.
namespace rhgpu {

constexpr int CoefficientWidth = 96; // 3*N-3 <= 93 for the curvature denominator.
constexpr int CoefficientHeight = 16;
constexpr int CurvatureMetaColumn = 95;
constexpr int CurvatureMetaRow = 15;

// u=y/4, P(u)=T(4u). The power-of-two variable change keeps the high-degree
// coefficients representable in float32 without deleting small terms.
// H=P*P'+u*(P*P''-P'^2), J=H'*P-2*H*P'. Consequently
// F(y)/F'(y)=4*H(u)*P(u)/J(u) for F=(y*T'/T)'.
// P/H/J are normalized separately by exact powers of two. Their quotient
// scale is included explicitly; no denominator or coefficient is clipped.
struct StableCurvatureData {
    int degree=0;
    long double variableScale=4.L,thetaScale=1.L,hScale=1.L,jScale=1.L,quotientScale=1.L;
    std::vector<long double> theta,h,j;
    bool valid=false;
    std::string reason;
};

inline StableCurvatureData makeStableCurvatureData(int degree) {
    StableCurvatureData data;
    data.degree=degree;
    if(degree<4||degree>32){data.reason="Unsupported stable curvature degree";return data;}
    data.theta.resize(degree+1);
    for(int k=0;k<=degree;++k)
        data.theta[k]=std::ldexp(static_cast<long double>(ff_rh::THETA_COEFFICIENTS[k]),2*k);
    data.h.assign(2*degree-1,0.L);
    // Symmetrizing TT'+u(TT''-T'^2) gives positive pair coefficients:
    // h[k]=sum_{i<j,i+j=k+1}(i-j)^2*a[i]*a[j].
    for(int i=0;i<=degree;++i)
        for(int j=i+1;j<=degree;++j) {
            const long double difference=j-i;
            data.h[i+j-1]+=difference*difference*data.theta[i]*data.theta[j];
        }
    data.j.assign(3*degree-2,0.L);
    std::vector<long double> compensation(data.j.size(),0.L),bounds(data.j.size(),0.L);
    for(std::size_t p=0;p<data.h.size();++p)
        for(int k=0;k<=degree;++k) {
            if(p==0&&k==0)continue;
            const std::size_t index=p+std::size_t(k)-1;
            const long double term=(static_cast<long double>(p)-2.L*k)*data.h[p]*data.theta[k];
            const long double sum=data.j[index],next=sum+term;
            compensation[index]+=std::abs(sum)>=std::abs(term)?(sum-next)+term:(term-next)+sum;
            data.j[index]=next;bounds[index]+=std::abs(term);
        }
    for(std::size_t k=0;k<data.j.size();++k) {
        data.j[k]+=compensation[k];
        // A badly cancelled host coefficient cannot be trusted after packing.
        // Decline only this optional payload; baseline GPU coefficients remain.
        if(bounds[k]>0.L&&std::abs(data.j[k])<=128.L*std::numeric_limits<long double>::epsilon()*bounds[k]) {
            data.reason="Stable curvature coefficient construction is ill conditioned";return data;
        }
    }
    auto normalize=[&](std::vector<long double>& polynomial,long double& scale) {
        long double largest=0.L;
        for(long double coefficient:polynomial) {
            if(!std::isfinite(coefficient))return false;
            largest=std::max(largest,std::abs(coefficient));
        }
        if(largest==0.L)return false;
        int exponent=0;std::frexp(largest,&exponent);scale=std::ldexp(1.L,exponent);
        for(auto& coefficient:polynomial) {
            coefficient/=scale;
            if(coefficient!=0.L&&std::abs(coefficient)<static_cast<long double>(std::numeric_limits<float>::min()))
                return false;
            const float packed=static_cast<float>(coefficient);
            if(!std::isfinite(packed)||(coefficient!=0.L&&packed==0.f))return false;
        }
        return true;
    };
    if(!normalize(data.theta,data.thetaScale)||!normalize(data.h,data.hScale)||!normalize(data.j,data.jScale)) {
        data.reason="Stable curvature coefficients exceed normal float32 range";return data;
    }
    data.quotientScale=data.variableScale*data.hScale*data.thetaScale/data.jScale;
    const float packedScale=static_cast<float>(data.quotientScale);
    if(!std::isfinite(packedScale)||packedScale==0.f||std::abs(packedScale)<std::numeric_limits<float>::min()) {
        data.reason="Stable curvature quotient scale exceeds normal float32 range";return data;
    }
    data.valid=true;
    return data;
}

struct alignas(16) Constants {
    std::uint32_t dispatch[4], selection[4], modes[4], style[4];
    float view[4], parameters[4], model[4], dynamics[4], paint[4], extra[4];
};
static_assert(sizeof(Constants) == 160, "RH compute constants must match HLSL packing");

inline Constants makeConstants(const amp::Scene& scene, int width, int height,
                               double aspect, int samples) {
    const amp::State& s = scene.state;
    const auto& e = s.rhExperiment;
    const auto& camera = s.cameras[std::clamp(s.plane, 0, 2)];
    Constants c{};
    c.dispatch[0] = static_cast<std::uint32_t>(width);
    c.dispatch[1] = static_cast<std::uint32_t>(height);
    c.dispatch[2] = 0;
    c.dispatch[3] = static_cast<std::uint32_t>(height);
    c.selection[0] = static_cast<std::uint32_t>(s.plane);
    c.selection[1] = static_cast<std::uint32_t>(s.generation - 1);
    c.selection[2] = static_cast<std::uint32_t>(8 << (std::clamp(s.spectralOrder, 1, 3) - 1));
    c.selection[3] = static_cast<std::uint32_t>(s.iterations);
    c.modes[0] = static_cast<std::uint32_t>(e.map);
    c.modes[1] = static_cast<std::uint32_t>(e.chart);
    c.modes[2] = static_cast<std::uint32_t>(e.coloring);
    c.modes[3] = static_cast<std::uint32_t>(e.power);
    c.style[0] = static_cast<std::uint32_t>(s.palette);
    c.style[1] = static_cast<std::uint32_t>(samples);
    c.style[2] = 0; // Engine enables metric writes only for eligible cache misses.
    c.view[0] = static_cast<float>(camera.x);
    c.view[1] = static_cast<float>(camera.y);
    c.view[2] = static_cast<float>(camera.span);
    c.view[3] = static_cast<float>(aspect);
    c.parameters[0] = static_cast<float>(s.c.real());
    c.parameters[1] = static_cast<float>(s.c.imag());
    c.parameters[2] = static_cast<float>(s.alpha.real());
    c.parameters[3] = static_cast<float>(s.alpha.imag());
    c.model[0] = static_cast<float>(s.seed.imag());
    c.model[1] = static_cast<float>(s.spectralD);
    c.model[2] = static_cast<float>(scene.masterModel ? scene.masterModel->level : s.rhLevel);
    c.model[3] = static_cast<float>(e.twist);
    c.dynamics[0] = static_cast<float>(e.gain);
    c.dynamics[1] = static_cast<float>(e.relaxation);
    c.dynamics[2] = static_cast<float>(e.memory);
    c.dynamics[3] = static_cast<float>(e.trapRadius);
    c.paint[0] = static_cast<float>(e.trapAngle);
    c.paint[1] = static_cast<float>(e.exposure);
    c.paint[2] = static_cast<float>(e.frequency);
    c.paint[3] = static_cast<float>(s.hue);
    c.extra[0] = static_cast<float>(s.bands);
    c.extra[1] = 0.f; // xStart; the host narrows this for each dispatch tile.
    c.extra[2] = static_cast<float>(width); // xCount; full width by default.
    c.extra[3] = static_cast<float>(e.blend);
    return c;
}

struct Data {
    std::vector<std::array<float, 4>> texels;
    bool valid = false;
    std::string reason;
};

inline Data makeData(const amp::State& s) {
    Data data;
    data.texels.resize(CoefficientWidth * CoefficientHeight, std::array<float, 4>{});
    if (s.family != 20 || s.generation < 1 || s.generation > rhfresh::CandidateCount) {
        data.reason = "The RH GPU renderer requires a registered RH candidate.";
        return data;
    }
    if (s.spectralOrder < 1 || s.spectralOrder > 3 || !std::isfinite(s.seed.real())
        || s.seed.real() <= 0.) {
        data.reason = "Invalid finite Theta degree or radius for the RH GPU renderer.";
        return data;
    }
    const int degree = 8 << (s.spectralOrder - 1);
    // Shader Model 5 hardware may flush subnormals to zero. Reject them as
    // well as actual cast underflow rather than silently dropping a term.
    auto put = [&](int row, int column, long double value) {
        if (!std::isfinite(value)
            || std::abs(value) > static_cast<long double>(std::numeric_limits<float>::max())) {
            data.reason = "RH coefficient exceeds the finite float32 range; using CPU.";
            return false;
        }
        if (value != 0.L
            && std::abs(value) < static_cast<long double>(std::numeric_limits<float>::min())) {
            data.reason = "RH coefficient would underflow float32 shader arithmetic; using CPU.";
            return false;
        }
        const float converted = static_cast<float>(value);
        if (!std::isfinite(converted) || (value != 0.L && converted == 0.f)) {
            data.reason = "RH coefficient cannot be represented in float32; using CPU.";
            return false;
        }
        data.texels[static_cast<std::size_t>(row) * CoefficientWidth + column][0] = converted;
        return true;
    };
    for (int n = 0; n <= degree; ++n)
        if (!put(0, n, ff_rh::THETA_COEFFICIENTS[n])) return data;

    const int candidate = s.generation - 1;
    // Unused texture cell (row 15, col 94) transports R to the shader
    // without altering the stable 160-byte RH constant-buffer ABI.
    data.texels[std::size_t(15)*CoefficientWidth+94][0]=static_cast<float>(s.seed.real());
    if(candidate==0) {
        const auto stable=makeStableCurvatureData(degree);
        if(stable.valid) {
            for(std::size_t k=0;k<stable.theta.size();++k)data.texels[k][1]=static_cast<float>(stable.theta[k]);
            for(std::size_t k=0;k<stable.h.size();++k)data.texels[k][2]=static_cast<float>(stable.h[k]);
            for(std::size_t k=0;k<stable.j.size();++k)data.texels[k][3]=static_cast<float>(stable.j[k]);
            data.texels[std::size_t(CurvatureMetaRow)*CoefficientWidth+CurvatureMetaColumn]=
                {{1.f,static_cast<float>(stable.variableScale),static_cast<float>(stable.quotientScale),0.f}};
        }
        // A zero readiness marker delegates the optional route to CPU. It
        // never invalidates the unchanged reference shader's coefficient data.
    }
    const bool angular = candidate == 1 || candidate == 2 || candidate == 4
        || candidate == 8 || candidate == 9 || candidate == 10 || candidate == 12;
    if (angular) {
        for (int radial = 0; radial <= 2; ++radial) {
            auto coefficients = rhfresh::detail::angularAmplitudes(degree, s.seed.real(), radial);
            for (int derivative = 0; derivative <= 4; ++derivative) {
                const int row = 1 + radial * 5 + derivative;
                for (std::size_t n = 0; n < coefficients.size(); ++n)
                    if (!put(row, static_cast<int>(n), coefficients[n])) return data;
                if (derivative != 4)
                    coefficients = rhfresh::detail::chebyshevDerivative(coefficients);
            }
        }
    }
    data.valid = true;
    data.reason = "Finite Theta coefficients prepared for float32 GPU preview.";
    return data;
}

} // namespace rhgpu
