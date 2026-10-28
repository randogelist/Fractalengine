#pragma once
// Direct3D 11 / Shader Model 5 implementation of the RH Lab finite models.
// Float arithmetic is guarded, never repaired by changing a denominator.
// 0xffffffff asks the host to recompute this entire pixel with the CPU model.
namespace rhgpu {
inline constexpr const char rhLabShader[] = R"RHLABHLSL(
cbuffer LabConstants : register(b0) {
    uint4 dispatchCfg;   // width, height, first row, row count
    uint4 selectionCfg;  // plane, candidate (zero based), degree, iteration budget
    uint4 modeCfg;       // map, chart, colouring, power
    uint4 styleCfg;      // palette, samples per axis, write metrics, stable-curvature bit0
    float4 viewCfg;      // camera x, camera y, horizontal span, height/width
    float4 parameterCfg; // c.real, c.imag, lambda.real, lambda.imag
    float4 modelCfg;     // scale b, shift h, literal level, Cayley twist
    float4 dynamicsCfg;  // gain, relaxation, memory, circle radius
    float4 paintCfg;     // cross angle, exposure, frequency, hue
    float4 extraCfg;     // bands, first column, column count, carrier blend
};
Texture2D<float4> coefficientTable : register(t0);
Texture2D<float4> cachedMetrics : register(t1);
RWTexture2D<uint> outputPixels : register(u0);
RWTexture2D<float4> outputMetrics : register(u1);
static const float PI = 3.14159265358979323846;
static const float TAU = 6.28318530717958647692;
static const float CONDITION_LIMIT = 0.00003;
static const uint REPAIR_PIXEL = 0xffffffffu;

bool cfinite(float2 z) { return all(isfinite(z)); }
bool czero(float2 z) { return all(z == float2(0.,0.)); }
float carg(float2 z) { return czero(z) ? 0. : atan2(z.y,z.x); }
float2 cmul(float2 a,float2 b) {
    return float2(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);
}
float cabs(float2 z) {
    float largest = max(abs(z.x),abs(z.y));
    if(largest == 0.) return 0.;
    float2 scaled = z/largest;
    return largest*sqrt(dot(scaled,scaled));
}
float logRadius(float2 z) {
    float largest = max(abs(z.x),abs(z.y));
    if(largest == 0.) return -3.402823466e38;
    float2 scaled = z/largest;
    return log(largest)+.5*log(dot(scaled,scaled));
}
float logOnePlusRadius(float2 z) {
    if(czero(z)) return 0.;
    float lr = logRadius(z);
    return lr > 36. ? lr : log(1.+exp(lr));
}
// Smith's complex division avoids squaring a very small or large denominator.
// denominatorBound is the sum of magnitudes of terms that formed it. A small
// denominator caused by cancellation is delegated to the precise CPU path.
float2 cdivide(float2 numerator,float2 denominator,float denominatorBound,inout bool good) {
    if(!good || !cfinite(numerator) || !cfinite(denominator) || !isfinite(denominatorBound)) {
        good=false; return float2(0.,0.);
    }
    if(czero(denominator)) {
        // Float cancellation cannot certify a true common zero. The CPU
        // applies the model's exact 0/0 continuation or singularity status.
        good=false; return float2(0.,0.);
    }
    float magnitude=cabs(denominator);
    if(magnitude < CONDITION_LIMIT*denominatorBound || magnitude < 1.e-30) {
        good=false; return float2(0.,0.);
    }
    float2 value;
    if(abs(denominator.x)>=abs(denominator.y)) {
        float ratio=denominator.y/denominator.x;
        float normal=denominator.x+denominator.y*ratio;
        if(!isfinite(normal)||normal==0.){good=false;return float2(0.,0.);}
        value=float2((numerator.x+numerator.y*ratio)/normal,
                     (numerator.y-numerator.x*ratio)/normal);
    } else {
        float ratio=denominator.x/denominator.y;
        float normal=denominator.y+denominator.x*ratio;
        if(!isfinite(normal)||normal==0.){good=false;return float2(0.,0.);}
        value=float2((numerator.x*ratio+numerator.y)/normal,
                     (numerator.y*ratio-numerator.x)/normal);
    }
    if(!cfinite(value))good=false;
    return value;
}
float coefficient(int index,int rowIndex) {
    return coefficientTable.Load(int3(index,rowIndex,0)).x;
}
float2 cpow(float2 z,uint exponent) {
    float2 value=cmul(z,z);
    [loop] for(uint k=2;k<exponent;++k)value=cmul(value,z);
    return value;
}

// Exactly the CPU polynomial jets: regular Horner inside the unit disk;
// reversed Horner outside, returning T^(k)(y)/y^(degree-k).
void polynomialJets(float2 y,bool reversed,int maxOrder,out float2 jets[5],inout bool good) {
    float2 argument=y;
    if(reversed)argument=cdivide(float2(1.,0.),y,cabs(y),good);
    [unroll] for(int clear=0;clear<5;++clear)jets[clear]=float2(0.,0.);
    [loop] for(int order=0;order<=maxOrder;++order) {
        float2 value=float2(0.,0.);
        float absoluteBound=0.;
        float argumentSize=cabs(argument);
        if(reversed) {
            [loop] for(int j=order;j<=(int)selectionCfg.z;++j) {
                float a=coefficient(j,0);
                [loop] for(int factor=0;factor<order;++factor)a*=j-factor;
                value=cmul(value,argument)+float2(a,0.);
                absoluteBound=absoluteBound*argumentSize+abs(a);
            }
        } else {
            [loop] for(int j=(int)selectionCfg.z;j>=order;--j) {
                float a=coefficient(j,0);
                [loop] for(int factor=0;factor<order;++factor)a*=j-factor;
                value=cmul(value,argument)+float2(a,0.);
                absoluteBound=absoluteBound*argumentSize+abs(a);
            }
        }
        jets[order]=value;
        if(!cfinite(value)||!isfinite(absoluteBound))good=false;
        // A theta zero gives correction zero in the exact model, but a
        // float-cancelled value requires recomputation to identify that zero.
        if(order==0&&absoluteBound>0.&&cabs(value)<CONDITION_LIMIT*absoluteBound)good=false;
    }
}
// Shared-normalization Horner for finite-shift and projective targets.
// Q(x)=T(M*x)/M^N and output jet[k]=T^(k)(w)/M^N.
void boundedThetaJets(float2 w,float M,int maximumOrder,out float2 jets[5],
                      out float bounds[5],inout bool good) {
    [unroll] for(int clear=0;clear<5;++clear){jets[clear]=float2(0.,0.);bounds[clear]=0.;}
    if(!cfinite(w)||!isfinite(M)||M<1.){good=false;return;}
    float2 x=w/M;
    float xSize=cabs(x);
    float inverse=1./M, inversePower=1.;
    [loop] for(int j=(int)selectionCfg.z;j>=0;--j) {
        [loop] for(int k=maximumOrder;k>=1;--k) {
            jets[k]=cmul(x,jets[k])+(float)k*jets[k-1];
            bounds[k]=xSize*bounds[k]+(float)k*bounds[k-1];
        }
        float term=coefficient(j,0)*inversePower;
        jets[0]=cmul(x,jets[0])+float2(term,0.);
        bounds[0]=xSize*bounds[0]+abs(term);
        inversePower*=inverse;
    }
    float derivativeScale=inverse;
    [loop] for(int k=1;k<=maximumOrder;++k) {
        jets[k]*=derivativeScale;
        bounds[k]*=derivativeScale;
        derivativeScale*=inverse;
    }
    [loop] for(int k=0;k<=maximumOrder;++k)
        if(!cfinite(jets[k])||!isfinite(bounds[k]))good=false;
    if(bounds[0]>0.&&cabs(jets[0])<CONDITION_LIMIT*bounds[0])good=false;
}
float2 angularValue(int radialOrder,int angularOrder,float2 y,inout bool good) {
    float2 newer=float2(0.,0.),older=float2(0.,0.);
    float newerBound=0.,olderBound=0.;
    float argumentSize=cabs(y);
    int rowIndex=1+radialOrder*5+angularOrder;
    // Rows are zero padded to 33 columns; after k derivatives the degree
    // is N-k, exactly as in the CPU's trimmed Chebyshev polynomial.
    int degree=max(0,(int)selectionCfg.z-angularOrder);
    [loop] for(int k=degree;k>=1;--k) {
        float a=coefficient(k,rowIndex);
        float2 current=2.*cmul(y,newer)-older+float2(a,0.);
        float currentBound=2.*argumentSize*newerBound+olderBound+abs(a);
        older=newer;newer=current;
        olderBound=newerBound;newerBound=currentBound;
    }
    float2 value=cmul(y,newer)-older+float2(coefficient(0,rowIndex),0.);
    float absoluteBound=argumentSize*newerBound+olderBound+abs(coefficient(0,rowIndex));
    if(!cfinite(value)||!isfinite(absoluteBound)
        ||(absoluteBound>0.&&cabs(value)<CONDITION_LIMIT*absoluteBound))good=false;
    return value;
}
// Optional curvature experiment. Coefficients were constructed in long
// double from the same binary64 finite model, in the variable u=y/4.
// Each row is power-of-two normalized; no term is dropped to fit float32.
// Both Horner branches keep |argument|<=1 to avoid high powers/overflow.
float2 curvaturePolynomial(int channel,int degree,float2 argument,bool reversed,
                           out float absoluteBound,inout bool good) {
    float2 value=float2(0.,0.);
    absoluteBound=0.;
    float argumentSize=cabs(argument);
    [loop] for(int n=0;n<=degree;++n) {
        int index=reversed?n:degree-n;
        float a=coefficientTable.Load(int3(index,0,0))[channel];
        value=cmul(value,argument)+float2(a,0.);
        absoluteBound=absoluteBound*argumentSize+abs(a);
    }
    if(!cfinite(value)||!isfinite(absoluteBound))good=false;
    return value;
}
float2 stableCurvatureCorrection(float2 y,inout bool good) {
    float4 meta=coefficientTable.Load(int3(95,15,0));
    if(meta.x!=1.||meta.y!=4.||!all(isfinite(meta))||meta.z<=0.) {
        good=false;return float2(0.,0.);
    }
    float2 u=y/meta.y;
    bool reversed=cabs(u)>1.;
    float2 argument=u;
    if(reversed)argument=cdivide(float2(1.,0.),u,cabs(u),good);
    int degree=(int)selectionCfg.z;
    float thetaBound,hBound,jBound;
    float2 theta=curvaturePolynomial(1,degree,argument,reversed,thetaBound,good);
    float2 h=curvaturePolynomial(2,2*degree-2,argument,reversed,hBound,good);
    float2 j=curvaturePolynomial(3,3*degree-3,argument,reversed,jBound,good);
    // Keep the original theta-root guard. Near a true or rounded common
    // zero, the CPU decides the continuation; J's own Horner bound guards
    // cancellation at a Newton pole without subtracting logarithmic jets.
    if(thetaBound>0.&&cabs(theta)<CONDITION_LIMIT*thetaBound)good=false;
    float2 value=cdivide(cmul(h,theta),j,jBound,good);
    // deg(H*P)=3N-2, deg(J)=3N-3. Reversed evaluation therefore
    // contributes exactly one u, followed by the explicit row/variable scale.
    if(reversed)value=cmul(u,value);
    value*=meta.z;
    if(!cfinite(value))good=false;
    return value;
}
// C(y)=degreeBase+offset and C' at the actual supplied argument. For
// |y|>1 use B(u)=u^N*T(1/u), u=1/y: C=N-u*B'/B. Keeping N separate
// avoids rounding away the small dilation gap in the difference C(2y)-C(y).
void logarithmicCountAndDerivative(float2 y,out float2 offset,out int degreeBase,out float2 derivative,
                                  out float derivativeBound,inout bool good) {
    bool reversed=cabs(y)>1.;
    degreeBase=0;
    if(reversed) {
        degreeBase=(int)selectionCfg.z;
        float2 inverse=cdivide(float2(1.,0.),y,cabs(y),good);
        float inverseSize=cabs(inverse);
        float2 inverseSquared=cmul(inverse,inverse);
        // Flushed/subnormal inverse powers must not create a false zero gap
        // or zero derivative. The CPU retains these values in wider arithmetic.
        if(inverseSize<1.175494351e-38||cabs(inverseSquared)<1.175494351e-38)good=false;
        float2 inverseJets[3];
        [loop] for(int order=0;order<=2;++order) {
            float2 term=float2(0.,0.);float bound=0.;
            [loop] for(int j=degreeBase;j>=order;--j) {
                float a=coefficient(degreeBase-j,0);
                [loop] for(int factor=0;factor<order;++factor)a*=j-factor;
                term=cmul(term,inverse)+float2(a,0.);
                bound=bound*inverseSize+abs(a);
            }
            inverseJets[order]=term;
            if(!cfinite(term)||!isfinite(bound))good=false;
            if(order==0&&bound>0.&&cabs(term)<CONDITION_LIMIT*bound)good=false;
        }
        float2 r=cdivide(inverseJets[1],inverseJets[0],cabs(inverseJets[0]),good);
        float2 s=cdivide(inverseJets[2],inverseJets[0],cabs(inverseJets[0]),good);
        float2 rr=cmul(r,r),logSecond=s-rr;
        offset=-cmul(inverse,r);
        derivative=cmul(inverseSquared,r+cmul(inverse,logSecond));
        derivativeBound=cabs(inverseSquared)*(cabs(r)+inverseSize*(cabs(s)+cabs(rr)));
        if(cabs(offset)<1.175494351e-38||cabs(derivative)<1.175494351e-38)good=false;
    } else {
        float2 jets[5];polynomialJets(y,false,2,jets,good);
        float2 a=cdivide(jets[1],jets[0],cabs(jets[0]),good);
        float2 b=cdivide(jets[2],jets[0],cabs(jets[0]),good);
        float2 aa=cmul(a,a),l2=b-aa;
        float l2Bound=cabs(b)+cabs(aa);
        offset=cmul(y,a);
        derivative=a+cmul(y,l2);
        derivativeBound=cabs(a)+cabs(y)*l2Bound;
    }
    if(!cfinite(offset)||!cfinite(derivative)||!isfinite(derivativeBound))good=false;
}
// The returned value is the NEGATIVE Newton correction, matching
// rhmaster::CandidateChannel::response, including division by model scale b.
float2 candidateResponse(float2 z,inout bool good) {
    float2 y=modelCfg.x*z+float2(modelCfg.y,0.);
    bool reversed=cabs(y)>1.;
    uint candidateIndex=selectionCfg.y;
    float target=modelCfg.z;
    float2 value=float2(0.,0.);
    if(candidateIndex==0&&(styleCfg.w&1u)!=0) {
        value=-stableCurvatureCorrection(y,good)/modelCfg.x;
        if(!cfinite(value))good=false;
        return value;
    }
    if(candidateIndex==0||candidateIndex==3||candidateIndex==6||candidateIndex==7) {
        // Curvature needs derivatives through order3, integer levels order2.
        int maxOrder=candidateIndex==3?2:(candidateIndex==0?3:4);
        float2 jets[5];polynomialJets(y,reversed,maxOrder,jets,good);
        if(!good)return value;
        if(czero(jets[0]))return value;
        float2 a=cdivide(jets[1],jets[0],cabs(jets[0]),good);
        float2 b=cdivide(jets[2],jets[0],cabs(jets[0]),good);
        float2 aa=cmul(a,a);
        float2 l2=b-aa;
        float l2Bound=cabs(b)+cabs(aa);
        float2 q=reversed?a+l2:a+cmul(y,l2);
        if(candidateIndex==3) {
            float qBound=reversed?cabs(a)+l2Bound:cabs(a)+cabs(y)*l2Bound;
            value=cdivide((reversed?a:cmul(y,a))-float2(target,0.),q,qBound,good);
        } else {
            float2 c=cdivide(jets[3],jets[0],cabs(jets[0]),good);
            float2 ab=cmul(a,b),aaa=cmul(aa,a);
            float2 l3=c-3.*ab+2.*aaa;
            float l3Bound=cabs(c)+3.*cabs(ab)+2.*cabs(aaa);
            float2 qp=reversed?2.*l2+l3:2.*l2+cmul(y,l3);
            float qpBound=reversed?2.*l2Bound+l3Bound:2.*l2Bound+cabs(y)*l3Bound;
            if(candidateIndex==0)value=cdivide(q,qp,qpBound,good);
            else {
                float2 d=cdivide(jets[4],jets[0],cabs(jets[0]),good);
                float2 l4=d-4.*cmul(a,c)-3.*cmul(b,b)+12.*cmul(aa,b)-6.*cmul(aa,aa);
                float l4Bound=cabs(d)+4.*cabs(cmul(a,c))+3.*cabs(cmul(b,b))
                    +12.*cabs(cmul(aa,b))+6.*cabs(cmul(aa,aa));
                float2 qpp=reversed?3.*l3+l4:3.*l3+cmul(y,l4);
                float qppBound=reversed?3.*l3Bound+l4Bound:3.*l3Bound+cabs(y)*l4Bound;
                if(candidateIndex==6)value=cdivide(qp,qpp,qppBound,good);
                else {
                    float2 numerator=q+(reversed?qp:cmul(y,qp));
                    float2 denominator=2.*qp+(reversed?qpp:cmul(y,qpp));
                    float denominatorBound=2.*qpBound+(reversed?1.:cabs(y))*qppBound;
                    value=cdivide(numerator,denominator,denominatorBound,good);
                }
            }
        }
        if(reversed)value=cmul(y,value);
    } else if(candidateIndex==5) {
        float2 v=-y;
        float2 jets[5];polynomialJets(v,reversed,3,jets,good);
        if(!good)return value;
        float2 t01=cmul(jets[0],jets[1]),t02=cmul(jets[0],jets[2]);
        float2 t11=cmul(jets[1],jets[1]),t03=cmul(jets[0],jets[3]);
        float2 t12=cmul(jets[1],jets[2]);
        float2 h,hp;
        float hBound,hpBound;
        if(reversed) {
            h=t01+t02-t11;hp=2.*t02+t03-t12;
            hBound=cabs(t01)+cabs(t02)+cabs(t11);
            hpBound=2.*cabs(t02)+cabs(t03)+cabs(t12);
            float2 inverse=cdivide(float2(1.,0.),v,cabs(v),good);
            float2 inversePower=float2(1.,0.);
            [loop] for(uint j=0;j<2*selectionCfg.z;++j)inversePower=cmul(inversePower,inverse);
            value=cmul(v,cdivide(-2.*h-target*inversePower,2.*(h+hp),2.*(hBound+hpBound),good));
        } else {
            h=t01+cmul(v,t02-t11);hp=2.*t02+cmul(v,t03-t12);
            hBound=cabs(t01)+cabs(v)*(cabs(t02)+cabs(t11));
            hpBound=2.*cabs(t02)+cabs(v)*(cabs(t03)+cabs(t12));
            value=cdivide(2.*cmul(y,h)-float2(target,0.),2.*h-2.*cmul(y,hp),
                2.*hBound+2.*cabs(y)*hpBound,good);
        }
    } else if(candidateIndex==1) {
        float2 u=angularValue(0,1,y,good),ux=angularValue(0,2,y,good);
        value=cdivide(u-float2(target,0.),ux,cabs(ux),good);
    } else if(candidateIndex==2) {
        float2 potential=angularValue(0,0,y,good),u=angularValue(0,1,y,good),ux=angularValue(0,2,y,good);
        float2 first=cmul(ux,potential),second=cmul(u,u);
        value=cdivide(cmul(u-target*potential,potential),first-second,cabs(first)+cabs(second),good);
    } else if(candidateIndex==4) {
        float2 ux=angularValue(0,2,y,good),uxx=angularValue(0,3,y,good);
        value=cdivide(ux,uxx,cabs(uxx),good);
    } else if(candidateIndex==8) {
        float2 ur=angularValue(1,1,y,good),urx=angularValue(1,2,y,good);
        value=cdivide(ur,urx,cabs(urx),good);
    } else if(candidateIndex==9) {
        float2 ux=angularValue(0,2,y,good),uxx=angularValue(0,3,y,good),uxxx=angularValue(0,4,y,good);
        float2 ur=angularValue(1,1,y,good),urx=angularValue(1,2,y,good),urxx=angularValue(1,3,y,good);
        float2 urr=angularValue(2,1,y,good),urrx=angularValue(2,2,y,good);
        float2 f1=cmul(urr,cmul(ux,ux));
        float2 f2=2.*cmul(urx,cmul(ur,ux));
        float2 f3=cmul(uxx,cmul(ur,ur));
        float2 d1=cmul(urrx,cmul(ux,ux));
        float2 d2=2.*cmul(urr,cmul(ux,uxx));
        float2 d3=2.*cmul(urxx,cmul(ur,ux));
        float2 d4=2.*cmul(cmul(urx,urx),ux);
        float2 d5=cmul(uxxx,cmul(ur,ur));
        value=cdivide(f1-f2+f3,d1+d2-d3-d4+d5,
            cabs(d1)+cabs(d2)+cabs(d3)+cabs(d4)+cabs(d5),good);
    } else if(candidateIndex==10) {
        // F=(E*E''-E'^2)/E^2-level, differentiated before division.
        float2 potential=angularValue(0,0,y,good),u=angularValue(0,1,y,good);
        float2 ux=angularValue(0,2,y,good),uxx=angularValue(0,3,y,good);
        float2 dFirst=cmul(potential,ux),dSecond=cmul(u,u);
        float2 dpFirst=cmul(potential,uxx),dpSecond=cmul(u,ux);
        float2 d=dFirst-dSecond,dp=dpFirst-dpSecond;
        float2 denominator=cmul(dp,potential)-2.*cmul(d,u);
        float denominatorBound=(cabs(dpFirst)+cabs(dpSecond))*cabs(potential)
            +2.*(cabs(dFirst)+cabs(dSecond))*cabs(u);
        value=cdivide(cmul(d-target*cmul(potential,potential),potential),
            denominator,denominatorBound,good);
    } else if(candidateIndex==11) {
        float2 offset,derivative,doubledOffset,doubledDerivative;
        int degreeBase,doubledDegreeBase;
        float derivativeBound,doubledDerivativeBound;
        logarithmicCountAndDerivative(y,offset,degreeBase,derivative,derivativeBound,good);
        logarithmicCountAndDerivative(2.*y,doubledOffset,doubledDegreeBase,doubledDerivative,doubledDerivativeBound,good);
        value=cdivide(doubledOffset-offset+float2((float)(doubledDegreeBase-degreeBase)-target,0.),
            2.*doubledDerivative-derivative,
            2.*doubledDerivativeBound+derivativeBound,good);
    } else if(candidateIndex==12) {
        // U=E'; F=U_R/U_x-level. All derivatives use the same finite lift.
        float2 ur=angularValue(1,1,y,good),ux=angularValue(0,2,y,good);
        float2 urx=angularValue(1,2,y,good),uxx=angularValue(0,3,y,good);
        float2 first=cmul(urx,ux),second=cmul(ur,uxx);
        value=cdivide(cmul(ur-target*ux,ux),first-second,cabs(first)+cabs(second),good);
    } else if(candidateIndex>=13&&candidateIndex<=17) {
        float2 a=y,b=y,c=y;
        int order=1;
        if(candidateIndex==13) {
            // Meta texel is reserved for the finite echo shift R.
            float echoRadius=coefficientTable.Load(int3(94,15,0)).x;
            if(!(echoRadius>0.)||!isfinite(echoRadius))good=false;
            a=y+float2(echoRadius,0.);b=y-float2(echoRadius,0.);
        } else if(candidateIndex==14)order=3;
        else if(candidateIndex==15){a=2.*y;b=.5*y;}
        else if(candidateIndex==16){b=-y;order=2;}
        else if(candidateIndex==17)order=4;
        float M=max(1.,cabs(a));
        if(candidateIndex==13||candidateIndex==15||candidateIndex==16)M=max(M,cabs(b));
        if(candidateIndex==13||candidateIndex==15)M=max(M,cabs(c));
        float2 A[5],B[5],C[5];float ab[5],bb[5],cb[5];
        boundedThetaJets(a,M,order,A,ab,good);
        if(!good)return value;
        float2 F=float2(0.,0.),Fp=float2(0.,0.);
        float numeratorBound=0.,denominatorBound=0.;
        if(candidateIndex==14) {
            float2 first=cmul(A[0],A[2]),second=cmul(A[1],A[1]);
            float2 derivFirst=cmul(A[0],A[3]),derivSecond=cmul(A[1],A[2]);
            F=first-second;Fp=derivFirst-derivSecond;
            numeratorBound=ab[0]*ab[2]+ab[1]*ab[1];
            denominatorBound=ab[0]*ab[3]+ab[1]*ab[2];
        } else if(candidateIndex==17) {
            float2 first=2.*cmul(A[1],A[3]),second=3.*cmul(A[2],A[2]);
            float2 derivFirst=2.*cmul(A[1],A[4]),derivSecond=4.*cmul(A[2],A[3]);
            F=first-second;Fp=derivFirst-derivSecond;
            numeratorBound=2.*ab[1]*ab[3]+3.*ab[2]*ab[2];
            denominatorBound=2.*ab[1]*ab[4]+4.*ab[2]*ab[3];
        } else {
            boundedThetaJets(b,M,order,B,bb,good);
            if(!good)return value;
            if(candidateIndex==16) {
                float2 first=cmul(A[1],B[0]),second=cmul(A[0],B[1]);
                float2 derivFirst=cmul(A[2],B[0]);
                float2 derivSecond=2.*cmul(A[1],B[1]);
                float2 derivThird=cmul(A[0],B[2]);
                F=first-second;Fp=derivFirst-derivSecond+derivThird;
                numeratorBound=ab[1]*bb[0]+ab[0]*bb[1];
                denominatorBound=ab[2]*bb[0]+2.*ab[1]*bb[1]+ab[0]*bb[2];
            } else {
                boundedThetaJets(c,M,order,C,cb,good);
                if(!good)return value;
                float firstRate=candidateIndex==13?1.:2.;
                float secondRate=candidateIndex==13?1.:.5;
                float2 first=cmul(A[0],B[0]),second=cmul(C[0],C[0]);
                float2 derivFirst=firstRate*cmul(A[1],B[0]);
                float2 derivSecond=secondRate*cmul(A[0],B[1]);
                float2 derivThird=2.*cmul(C[0],C[1]);
                F=first-second;Fp=derivFirst+derivSecond-derivThird;
                numeratorBound=ab[0]*bb[0]+cb[0]*cb[0];
                denominatorBound=firstRate*ab[1]*bb[0]+secondRate*ab[0]*bb[1]+2.*cb[0]*cb[1];
            }
        }
        if(!cfinite(F)||!cfinite(Fp)||!isfinite(numeratorBound)||!isfinite(denominatorBound))good=false;
        // Cancellation at a true target root is assigned to CPU for the
        // final precise classification. Never insert a fake denominator.
        if(numeratorBound>0.&&cabs(F)<CONDITION_LIMIT*numeratorBound)good=false;
        value=cdivide(F,Fp,denominatorBound,good);
    } else good=false;
    value=-value/modelCfg.x;
    if(!cfinite(value))good=false;
    return value;
}
float2 chartForward(float2 w,inout bool good) {
    if(modeCfg.y==0)return w;
    if(modeCfg.y==2) {
        float radius=exp(w.x);
        // A normal finite display coordinate can underflow in float exp;
        // let the CPU distinguish its positive value from a true zero.
        if(radius<1.175494351e-38){good=false;return float2(0.,0.);}
        float2 value=radius*float2(cos(w.y),sin(w.y));
        if(!cfinite(value))good=false;
        return value;
    }
    float2 reciprocal=cdivide(float2(1.,0.),w,cabs(w),good);
    if(modeCfg.y==1)return reciprocal;
    if(modeCfg.y==3)return .5*(w+reciprocal);
    good=false;return float2(0.,0.);
}

struct OrbitResult {
    float2 endpoint;
    uint iterations;
    uint escaped;
    uint detectedPeriod;
    float smooth;
    float structure;
    float phase;
    float periodResidual;
    float circleTrap;
    float crossTrap;
    float winding;
    float feedbackPhase;
    float feedbackMean;
    float meanLogRadius;
};
OrbitResult orbit(float2 initial,float2 parameter,float2 strength,inout bool good) {
    OrbitResult r;
    r.endpoint=initial;r.iterations=0;r.escaped=0;r.detectedPeriod=0;
    r.smooth=0.;r.structure=0.;r.phase=0.;r.periodResidual=1.e30;
    r.circleTrap=1.e30;r.crossTrap=1.e30;r.winding=0.;
    r.feedbackPhase=0.;r.feedbackMean=0.;r.meanLogRadius=0.;
    float2 z=initial,previous=initial,feedbackDirection=float2(0.,0.);
    float2 rotation=float2(cos(paintCfg.x),-sin(paintCfg.x));
    uint feedbackSamples=0,orbitSamples=0;
    float minimumTrap=1.e30;
    bool retainTail=modeCfg.z==6||(selectionCfg.x==2&&modeCfg.z==0);
    float2 tail[80]; // Only assigned/read when a finite-tail colouring needs it.
    [loop] for(uint n=0;n<=selectionCfg.w;++n) {
        if(!good||!cfinite(z)){good=false;break;}
        r.endpoint=z;r.iterations=n;r.smooth=(float)n;
        if(retainTail)tail[n%80]=z;
        if(modeCfg.z==0)minimumTrap=min(minimumTrap,cabs(z-float2(-.21,.37)));
        if(modeCfg.z==2)r.circleTrap=min(r.circleTrap,abs(cabs(z)-dynamicsCfg.w));
        if(modeCfg.z==3) {
            float2 rotated=cmul(z,rotation);
            r.crossTrap=min(r.crossTrap,min(abs(rotated.x),abs(rotated.y)));
        }
        if(modeCfg.z==4&&n>0&&!czero(z)&&!czero(previous)) {
            float difference=carg(z)-carg(previous);
            r.winding+=difference-TAU*round(difference/TAU);
        }
        if(modeCfg.z==7){r.meanLogRadius+=logOnePlusRadius(z);++orbitSamples;}
        float largest=max(abs(z.x),abs(z.y));
        if(largest>16.||dot(z,z)>256.) {
            r.escaped=1;
            r.smooth=max(0.,(float)n+1.-log(logRadius(z)/log(16.))/log(2.));
            break;
        }
        if(n==selectionCfg.w)break;
        float2 carrier=z;
        if(modeCfg.x==2||modeCfg.x==6||(modeCfg.x==7&&(n&1)))
            carrier=float2(z.x,-z.y);
        else if(modeCfg.x==10||modeCfg.x==11) {
            if(extraCfg.w==1.)carrier=float2(z.x,-z.y);
            else if(extraCfg.w!=0.)carrier=float2(z.x,(1.-2.*extraCfg.w)*z.y);
        }
        float2 correction=float2(0.,0.);
        if(!czero(strength)&&dynamicsCfg.x!=0.)correction=candidateResponse(carrier,good);
        if(!good)break;
        float2 feedback=strength.y==0.?strength.x*correction:cmul(strength,correction);
        if(dynamicsCfg.x!=1.)feedback*=dynamicsCfg.x;
        if(modeCfg.x==3&&(n&1))feedback=-feedback;
        if((modeCfg.x==9||modeCfg.x==11)&&modelCfg.w!=0.)
            feedback=cmul(feedback,float2(cos(modelCfg.w),sin(modelCfg.w)));
        if(modeCfg.z==5) {
            r.feedbackMean+=logOnePlusRadius(feedback);
            if(!czero(feedback)) {
                // The unit direction equals exp(i*arg(feedback)). Scaling
                // first keeps its norm finite for every finite feedback.
                float largestFeedback=max(abs(feedback.x),abs(feedback.y));
                float2 scaledFeedback=feedback/largestFeedback;
                feedbackDirection+=scaledFeedback/sqrt(dot(scaledFeedback,scaledFeedback));
            }
            ++feedbackSamples;
        }
        float2 target;
        if(modeCfg.x==1)target=z+parameter+feedback;
        else {
            float2 host=cpow(carrier,modeCfg.w)+parameter;
            if((modeCfg.x>=5&&modeCfg.x<=9)||modeCfg.x==11) {
                if(czero(feedback))target=host;
                else {
                    float2 denominator=float2(1.,0.)-feedback;
                    target=cmul(host,cdivide(float2(1.,0.)+feedback,denominator,1.+cabs(feedback),good));
                }
                if(modeCfg.x==8&&dynamicsCfg.z!=0.)target+=dynamicsCfg.z*(z-previous);
            } else {
                target=host+feedback;
                if(modeCfg.x==4)target+=dynamicsCfg.z*(z-previous);
            }
        }
        float2 next=dynamicsCfg.y==1.?target:z+dynamicsCfg.y*(target-z);
        if(!cfinite(next)){good=false;break;}
        previous=z;z=next;
    }
    if(modeCfg.z==0)r.structure=exp(-2.5*minimumTrap);
    r.phase=frac(.5+carg(r.endpoint)/TAU);
    if(feedbackSamples>0) {
        r.feedbackMean/=(float)feedbackSamples;
        r.feedbackPhase=frac(.5+carg(feedbackDirection)/TAU);
    }
    if(orbitSamples>0)r.meanLogRadius/=(float)orbitSamples;
    if(good&&retainTail&&!r.escaped&&r.iterations>=32) {
        uint available=min(r.iterations+1,80u);
        int latest=(int)(r.iterations%80);
        [loop] for(uint period=1;period<=16&&5*period<=available;++period) {
            if((modeCfg.x==3||modeCfg.x==7)&&(period&1))continue;
            float residual=0.;bool repeated=true;
            [loop] for(uint lag=0;lag<4*period;++lag) {
                float2 recent=tail[(latest+80-(int)lag)%80];
                float2 older=tail[(latest+80-(int)lag-(int)period)%80];
                float error=cabs(recent-older)/(1.+max(cabs(recent),cabs(older)));
                if(!isfinite(error)||error>1.e-7){repeated=false;break;}
                residual=max(residual,error);
            }
            if(repeated){r.detectedPeriod=period;r.periodResidual=residual;break;}
        }
    }
    return r;
}
float3 hsv(float hue,float saturation,float value) {
    float sector=frac(hue)*6.;
    int which=(int)sector;
    float fraction=sector-which;
    float p=value*(1.-saturation),q=value*(1.-saturation*fraction),t=value*(1.-saturation*(1.-fraction));
    if(which==0)return float3(value,t,p);
    if(which==1)return float3(q,value,p);
    if(which==2)return float3(p,value,t);
    if(which==3)return float3(p,q,value);
    if(which==4)return float3(t,p,value);
    return float3(value,p,q);
}
float3 originalColour(OrbitResult r,float bands) {
    bands=clamp(bands,0.,100.);
    float hue=paintCfg.w;
    if(!r.escaped)return hsv(hue+.61,.48,.019+.017*r.structure);
    float stripes=.5+.5*cos(bands*(.57*r.smooth+2.4*r.structure));
    float colourHue=styleCfg.x==1?hue+.66-.023*r.smooth
        :styleCfg.x==2?hue+.65*r.structure+.18*r.phase
        :hue+.57+.30*r.phase+.012*r.smooth;
    return hsv(colourHue,.74,.26+.72*stripes);
}
float3 lambdaDiagnosticColour(OrbitResult r,float bands) {
    bands=clamp(bands,0.,100.);
    if(r.escaped||r.detectedPeriod==0)return originalColour(r,bands);
    float precision=saturate((-log10(max(1.e-16,r.periodResidual))-7.)/7.);
    float periodHue=styleCfg.x==1?.11+.057*r.detectedPeriod:.52+.6180339887498948*r.detectedPeriod;
    float phaseHue=styleCfg.x==2?.09*r.phase:0.;
    float modulation=.94+.06*cos(bands*(TAU*r.phase+r.structure));
    float brightness=(.42+.30*precision+.08*r.structure)*modulation;
    return hsv(paintCfg.w+periodHue+phaseHue,.65,brightness);
}
float3 experimentColour(OrbitResult r) {
    float bands=clamp(extraCfg.x,0.,100.);
    float frequency=paintCfg.z;
    if(modeCfg.z==0) {
        float3 value=selectionCfg.x==2?lambdaDiagnosticColour(r,bands*frequency):originalColour(r,bands*frequency);
        return paintCfg.y==1.?value:saturate(value*paintCfg.y);
    }
    float metric=0.,phase=r.phase,saturation=.76,brightness=.7;
    if(modeCfg.z==1){metric=r.phase;phase=r.phase;brightness=r.escaped?.82:.72;}
    if(modeCfg.z==2) {
        metric=-log10(max(1.e-14,r.circleTrap));phase=.14*metric;
        brightness=.32+.62*exp(-frequency*r.circleTrap);
    }
    if(modeCfg.z==3) {
        metric=-log10(max(1.e-14,r.crossTrap));phase=.16*metric;
        brightness=.36+.60*exp(-frequency*r.crossTrap);
    }
    if(modeCfg.z==4){metric=r.winding/TAU;phase=.17*metric+.12*r.phase;brightness=.66;}
    if(modeCfg.z==5){metric=r.feedbackMean;phase=r.feedbackPhase+.13*metric;brightness=.38+.58*(1.-exp(-metric));}
    if(modeCfg.z==6) {
        if(r.detectedPeriod>0) {
            metric=(float)r.detectedPeriod;phase=.6180339887498948*metric;
            brightness=.60+.35*saturate((-log10(max(1.e-16,r.periodResidual))-7.)/7.);
        } else {
            metric=r.smooth;phase=.017*r.smooth+.1*r.phase;
            brightness=r.escaped?.62:.15;saturation=.48;
        }
    }
    if(modeCfg.z==7){metric=r.meanLogRadius;phase=.19*metric+.1*r.phase;brightness=.32+.63*(1.-exp(-metric));}
    float stripes=.5+.5*cos(bands*frequency*(.71*metric+TAU*r.phase));
    float value=saturate(paintCfg.y*brightness*(.66+.34*stripes));
    float offset=styleCfg.x==1?.07:styleCfg.x==2?.36:.59;
    return hsv(paintCfg.w+offset+frequency*phase,saturation,value);
}
// Cache contract: a single sample, one float4 per pixel, colours1..7.
// x is endpoint phase; y/z are the selected diagnostic; w packs period
// and escape as small exact integers. w<0 means CPU precision repair.
// Appearance fields never enter these statistics. Changed dynamics,
// geometry, traps or diagnostic invalidate the host's cache key.
float4 encodeMetrics(OrbitResult r) {
    float first=0.,second=0.;
    if(modeCfg.z==2)first=r.circleTrap;
    if(modeCfg.z==3)first=r.crossTrap;
    if(modeCfg.z==4)first=r.winding;
    if(modeCfg.z==5){first=r.feedbackMean;second=r.feedbackPhase;}
    if(modeCfg.z==6){first=r.periodResidual;second=r.smooth;}
    if(modeCfg.z==7)first=r.meanLogRadius;
    return float4(r.phase,first,second,(float)(r.detectedPeriod+32*r.escaped));
}
OrbitResult decodeMetrics(float4 m) {
    OrbitResult r;
    r.endpoint=float2(0.,0.);r.iterations=0;r.smooth=0.;
    r.escaped=(uint)(m.w/32.);r.detectedPeriod=(uint)m.w-32*r.escaped;
    r.structure=0.;r.phase=m.x;r.periodResidual=1.e30;
    r.circleTrap=1.e30;r.crossTrap=1.e30;r.winding=0.;
    r.feedbackPhase=0.;r.feedbackMean=0.;r.meanLogRadius=0.;
    if(modeCfg.z==2)r.circleTrap=m.y;
    if(modeCfg.z==3)r.crossTrap=m.y;
    if(modeCfg.z==4)r.winding=m.y;
    if(modeCfg.z==5){r.feedbackMean=m.y;r.feedbackPhase=m.z;}
    if(modeCfg.z==6){r.periodResidual=m.y;r.smooth=m.z;}
    if(modeCfg.z==7)r.meanLogRadius=m.y;
    return r;
}
uint packColour(float3 colour) {
    uint3 packed=(uint3)clamp(pow(max(colour,float3(0.,0.,0.)),.88)*255.+.5,0.,255.);
    return packed.z|(packed.y<<8)|(packed.x<<16);
}
float3 colourAt(float2 displayPosition,inout bool good,out float4 metrics) {
    float2 initial=float2(0.,0.);
    float2 parameter=parameterCfg.xy,strength=parameterCfg.zw;
    if(selectionCfg.x==0)parameter=displayPosition;
    if(selectionCfg.x==1)initial=chartForward(displayPosition,good);
    if(selectionCfg.x==2)strength=displayPosition;
    OrbitResult result=orbit(initial,parameter,strength,good);
    float3 value=experimentColour(result);
    metrics=float4(0.,0.,0.,-1.);
    if(styleCfg.z!=0&&good) {
        metrics=encodeMetrics(result);
        if(!all(isfinite(metrics)))good=false;
    }
    if(!all(isfinite(value)))good=false;
    return value;
}
[numthreads(8,4,1)]
void main(uint3 dispatchId : SV_DispatchThreadID) {
    uint x=dispatchId.x+(uint)extraCfg.y,y=dispatchCfg.z+dispatchId.y;
    if(dispatchId.x>=(uint)extraCfg.z||x>=dispatchCfg.x||y>=dispatchCfg.y||dispatchId.y>=dispatchCfg.w)return;
    bool good=true;
    float3 accumulated=float3(0.,0.,0.);
    float4 metrics=float4(0.,0.,0.,-1.);
    uint sampleCount=styleCfg.y;
    [loop] for(uint sy=0;sy<sampleCount;++sy) {
        [loop] for(uint sx=0;sx<sampleCount;++sx) {
            float2 subpixel=float2((float)x+((float)sx+.5)/(float)sampleCount,
                                   (float)y+((float)sy+.5)/(float)sampleCount);
            float2 position=float2(viewCfg.x+(subpixel.x/(float)dispatchCfg.x-.5)*viewCfg.z,
                viewCfg.y+(.5-subpixel.y/(float)dispatchCfg.y)*viewCfg.z*viewCfg.w);
            accumulated+=colourAt(position,good,metrics)/(float)(sampleCount*sampleCount);
            if(!good)break;
        }
        if(!good)break;
    }
    if(styleCfg.z!=0)outputMetrics[uint2(x,y)]=good?metrics:float4(0.,0.,0.,-1.);
    if(!good||sampleCount==0){outputPixels[uint2(x,y)]=REPAIR_PIXEL;return;}
    outputPixels[uint2(x,y)]=packColour(accumulated);
}
// Recolouring has no orbit evaluator or derivative work in its reachable
// call graph. Repair markers deliberately survive so the host recalculates
// those pixels with the current appearance settings on every cache hit.
[numthreads(8,4,1)]
void recolour(uint3 dispatchId : SV_DispatchThreadID) {
    uint x=dispatchId.x+(uint)extraCfg.y,y=dispatchCfg.z+dispatchId.y;
    if(dispatchId.x>=(uint)extraCfg.z||x>=dispatchCfg.x||y>=dispatchCfg.y||dispatchId.y>=dispatchCfg.w)return;
    float4 metrics=cachedMetrics.Load(int3(x,y,0));
    if(metrics.w<0.||!all(isfinite(metrics))){outputPixels[uint2(x,y)]=REPAIR_PIXEL;return;}
    float3 value=experimentColour(decodeMetrics(metrics));
    outputPixels[uint2(x,y)]=all(isfinite(value))?packColour(value):REPAIR_PIXEL;
}
)RHLABHLSL";
} // namespace rhgpu
