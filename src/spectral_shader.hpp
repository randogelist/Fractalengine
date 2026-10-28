#pragma once
// Source only. D3DCompile runs on the user's PC, lazily on family selection.
static const char spectralShader[] = R"SPECTRALHLSL(
cbuffer Settings : register(b0) {
    uint4 sizeRow;
    uint4 mode;
    float4 view;
    float4 param;
    float4 option;
    float4 refA;      // x = spectral half-spacing d; y = G2 middle-atom weight
    float4 refB;
    uint4 flags;
    uint4 garden;
    uint4 record;     // y = G2 spectral jet order (1..3)
};
RWTexture2D<uint> pixels : register(u0);

float2 cm(float2 a,float2 b) {
    return float2(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);
}
bool finitePair(float2 z) { return all(abs(z)<=3.0e38); }
float magnitude(float2 z) {
    float s=max(abs(z.x),abs(z.y));
    if(s==0)return 0;
    return s*sqrt(dot(z/s,z/s));
}
float2 divide(float2 a,float2 b) {
    float s=max(abs(b.x),abs(b.y));
    float2 q=b/s,u=a/s;
    return float2(dot(u,q),u.y*q.x-u.x*q.y)/dot(q,q);
}
float phase(float2 z) {
    if(all(z==float2(0,0)))return 0;
    return atan2(z.y,z.x);
}
float3 rainbow(float t) {
    return .52+.46*cos(6.28318530718*(float3(0,.33,.67)+t));
}
float escapeRadius(float2 c,float2 eta,float d) {
    return max(max(2.,2.*d),sqrt(2.*(magnitude(c)+5.*magnitude(eta)/9.)));
}

float3 spectralColourGen1(float2 p) {
    float2 z=mode.x==0?p:option.xy;
    float2 c=mode.x==1?p:param.xy;
    float2 eta=mode.x==2?p:param.zw;
    float d=refA.x,d2=d*d;
    bool quadratic=(d==0.)||all(eta==float2(0,0));
    float radius=escapeRadius(c,quadratic?float2(0,0):eta,d);
    // status: 0 finite iteration budget; 1 escaped; 2 numerical failure.
    uint status=0,steps=0;
    float trap=1.,smooth=0.;
    float2 phaseSum=0.;
    bool certifiedPole=false;
    if(!finitePair(z)||!finitePair(c)||!finitePair(eta)||!(d>=0.)||
       !(d2<=3.0e38)||!(radius<=3.0e38)||(!quadratic && d2==0.))
        status=2;
    [loop]for(uint i=0;i<=mode.w && status==0;i++) {
        float zm=magnitude(z);
        if(!(zm<=3.0e38)){status=2;break;}
        trap=min(trap,min(magnitude(z-float2(0,d)),magnitude(z+float2(0,d)))/(1.+zm));
        if(zm>0.)phaseSum+=z/zm;
        steps=i;
        if(zm>radius){status=1;break;}
        if(i==mode.w)break;
        float2 z2=cm(z,z),base=z2+c;
        if(!finitePair(base)){status=2;break;}
        if(quadratic){z=base;continue;}
        float2 den=z2+float2(d2,0),num=z2-float2(d2,0);
        float denScale=max(abs(den.x),abs(den.y));
        // Close to a pole, the reverse triangle inequality certifies that
        // the next iterate exceeds the escape radius, before any division.
        if(magnitude(den)<=1.0e-8*(1.+magnitude(z2)+d2)) {
            float numeratorSize=magnitude(eta)*d2*magnitude(num);
            float denominatorSize=dot(den,den);
            float bound=(radius+magnitude(base))*denominatorSize;
            // Small conservative margin reduces false certificates caused
            // by float rounding; this remains a numerical visualization.
            if(numeratorSize>bound*(1.+8.0e-6)) {
                status=1;steps=i+1;certifiedPole=true;break;
            }
            if(denScale==0.){status=2;break;}
        }
        float2 correction=cm(eta,divide(divide(num,den)*d2,den));
        float2 next=base+correction;
        if(!finitePair(next)){status=2;break;}
        z=next;
    }
    if(status==2)return float3(.35,.02,.20);
    if(status==1) {
        smooth=certifiedPole?float(steps):
            max(0.,float(steps)+1.-log(log(max(magnitude(z),2.)))/.69314718056);
    }
    float depth=-log(max(trap,1.e-12));
    float t=0.,light=0.;
    if(mode.z==0) {
        if(status!=1)return float3(0,0,0);
        t=option.z+.08+.075*option.w*smooth;
        light=.45+.50*(1.-exp(-.045*smooth));
    }else if(mode.z==1) {
        t=option.z+.18+.12*option.w*depth;
        light=.20+.70*(1.-exp(-.22*depth));
    }else {
        t=option.z+option.w*phase(phaseSum)/6.28318530718+
            .025*option.w*log(1.+float(steps));
        light=status==1?.40+.50*(1.-exp(-.06*smooth)):
                            .26+.38*(1.-exp(-.18*depth));
    }
    return rainbow(t)*light;
}

// G2: three spectral offsets (-d, .35d, d), weights ((1-m)/2,m,(1-m)/2).
// Actual positive atoms are q_i=1/2+t_i; x is the centered coordinate.
// The scaled polynomial form evaluates the removable atom singularities
// directly.  Its numerator and denominator squared have equal total degree.
bool jet2Data(float2 z,float d,float weight,uint order,
              out float2 numerator,out float2 denominator,out float majorant) {
    numerator=0.;denominator=0.;majorant=0.;
    float scale=max(magnitude(z),d);
    if(!(scale>0.)||!(scale<=3.0e38))return false;
    float beta=d/scale,beta2=beta*beta;
    if(!(beta>0.)||!(beta2>0.))return false;
    float2 u0=z/scale+float2(beta,0);
    float2 u1=z/scale-float2(.35*beta,0);
    float2 u2=z/scale-float2(beta,0);
    float2 v0=u0,v1=u1,v2=u2;
    if(order>=2){v0=cm(v0,u0);v1=cm(v1,u1);v2=cm(v2,u2);}
    if(order>=3){v0=cm(v0,u0);v1=cm(v1,u1);v2=cm(v2,u2);}
    float2 p0=cm(v0,u0),p1=cm(v1,u1),p2=cm(v2,u2);
    float outer=(1.-weight)*.5;
    float2 b0=outer*cm(p1,p2),b1=weight*cm(p0,p2),b2=outer*cm(p0,p1);
    denominator=b0+b1+b2;
    majorant=magnitude(b0)+magnitude(b1)+magnitude(b2);
    numerator=beta2*(
        outer*weight*1.8225*cm(cm(v0,v1),cm(p2,p2))+
        outer*outer*4.*cm(cm(v0,v2),cm(p1,p1))+
        weight*outer*.4225*cm(cm(v1,v2),cm(p0,p0)));
    return finitePair(numerator)&&finitePair(denominator)&&
           (majorant>0.)&&(majorant<=3.0e38);
}

float3 spectralColourGen2(float2 p) {
    float2 z=mode.x==0?p:option.xy;
    float2 c=mode.x==1?p:param.xy;
    float2 eta=mode.x==2?p:param.zw;
    float d=refA.x,weight=refA.y;
    uint order=record.y;
    bool quadratic=(d==0.)||all(eta==float2(0,0));
    // For s=1..3, |D_s(z)|<1 for |z|>=8d.  Outside this radius
    // |F(z)|>=|z|^2/2>|z|, so the orbit tends to infinity.
    float radius=quadratic?escapeRadius(c,float2(0,0),0.):
        max(max(2.,8.*d),sqrt(2.*(magnitude(c)+magnitude(eta))));
    uint status=0,steps=0;
    float trap=1.,smooth=0.;
    float2 phaseSum=0.;
    bool certifiedPole=false;
    if(!finitePair(z)||!finitePair(c)||!finitePair(eta)||
       !(d>=0.)||!(d<=.49)||!(weight>=.01)||!(weight<=.98)||
       order<1||order>3||!(radius<=3.0e38))status=2;
    [loop]for(uint i=0;i<=mode.w && status==0;i++) {
        float zm=magnitude(z);
        if(!(zm<=3.0e38)){status=2;break;}
        if(zm>0.)phaseSum+=z/zm;
        steps=i;
        if(zm>radius){status=1;break;}
        float2 numerator=0.,denominator=0.;
        float majorant=1.;
        if(!quadratic) {
            if(!jet2Data(z,d,weight,order,numerator,denominator,majorant)) {
                status=2;break;
            }
            // Cancellation ratio, not distance to a spectral atom.
            trap=min(trap,clamp(magnitude(denominator)/majorant,0.,1.));
        }
        if(i==mode.w)break;
        float2 base=cm(z,z)+c;
        if(!finitePair(base)){status=2;break;}
        if(quadratic){z=base;continue;}
        float denominatorSize=magnitude(denominator);
        float numeratorSize=magnitude(numerator);
        if(denominatorSize==0.) {
            if(numeratorSize>0.) {
                status=1;steps=i+1;certifiedPole=true;
            }else status=2; // An indeterminate 0/0 is never counted as escape.
            break;
        }
        // A small denominator is insufficient: use the reverse triangle
        // inequality with a float-rounding margin to certify the next escape.
        float left=magnitude(eta)*numeratorSize;
        float right=(radius+magnitude(base))*denominatorSize*denominatorSize;
        if(denominatorSize<=1.e-8*majorant && left>right*(1.+8.e-6)) {
            status=1;steps=i+1;certifiedPole=true;break;
        }
        float2 next=base+divide(cm(eta,divide(numerator,denominator)),denominator);
        if(!finitePair(next)) {
            // Evaluate this same escape test in log scale on overflow.
            float logTerm=log(magnitude(eta))+log(numeratorSize)-2.*log(denominatorSize);
            if(numeratorSize>0. && logTerm>log(radius+magnitude(base))+8.e-6) {
                status=1;steps=i+1;certifiedPole=true;
            }else status=2;
            break;
        }
        z=next;
    }
    if(status==2)return float3(.35,.02,.20);
    if(status==1) {
        smooth=certifiedPole?float(steps):
            max(0.,float(steps)+1.-log(log(max(magnitude(z),2.)))/.69314718056);
    }
    float depth=-log(max(trap,1.e-12));
    float t=0.,light=0.;
    if(mode.z==0) {
        if(status!=1)return float3(0,0,0);
        t=option.z+.08+.075*option.w*smooth;
        light=.45+.50*(1.-exp(-.045*smooth));
    }else if(mode.z==1) {
        t=option.z+.18+.12*option.w*depth;
        light=.20+.70*(1.-exp(-.22*depth));
    }else {
        t=option.z+option.w*phase(phaseSum)/6.28318530718+
            .025*option.w*log(1.+float(steps));
        light=status==1?.40+.50*(1.-exp(-.06*smooth)):
                            .26+.38*(1.-exp(-.18*depth));
    }
    return rainbow(t)*light;
}

float3 spectralColour(float2 p) {
    [branch]if(garden.y>=2)return spectralColourGen2(p);
    return spectralColourGen1(p);
}

[numthreads(16,8,1)]
void main(uint3 id:SV_DispatchThreadID) {
    uint x=id.x,y=id.y+sizeRow.z;
    if(x>=sizeRow.x||y>=sizeRow.y||id.y>=sizeRow.w)return;
    float3 col=0;
    [loop]for(uint sy=0;sy<flags.z;sy++)[loop]for(uint sx=0;sx<flags.z;sx++) {
        float2 uv=(float2(x,y)+(float2(sx,sy)+.5)/float(flags.z))/float2(sizeRow.xy);
        float2 p=view.xy+(uv-.5)*float2(view.z,-view.z*view.w);
        col+=spectralColour(p);
    }
    col/=float(flags.z*flags.z);
    uint3 rgb=(uint3)clamp(pow(max(col,0),.88)*255+.5,0,255);
    pixels[uint2(x,y)]=rgb.b|(rgb.g<<8)|(rgb.r<<16);
}
)SPECTRALHLSL";
