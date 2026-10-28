#pragma once
// Runtime-compiled Shader Model 5 compute shader. No external shader file is required.
static const char ampShader[] = R"AMPHLSL(
cbuffer Settings : register(b0) {
    uint4 sizeRow;       // width, height, first row, number of rows
    uint4 mode;          // plane, critical selection, palette, iteration budget
    float4 view;         // center real, center imag, horizontal span, height/width
    float4 param;        // c.real, c.imag, alpha.real, alpha.imag
    float4 option;       // seed.real, seed.imag, hue, band frequency
    float4 refA;
    float4 refB;
    uint4 flags;         // reference A valid, B valid, samples per axis, reserved
};
RWTexture2D<uint> pixels : register(u0);
float2 cm(float2 a,float2 b) { return float2(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x); }
float cd(float4 a,float4 b) { return length(cm(a.xy,b.zw)-cm(b.xy,a.zw)); }
float4 pt(float2 z) { return normalize(float4(z,1,0)); }
float4 advance(float4 q,float2 c,float2 alpha,out float ly) {
    if(dot(alpha,alpha)==0) {ly=log(1e-24);return pt(c);}
    float2 xx=cm(q.xy,q.xy),xy=cm(q.xy,q.zw),yy=cm(q.zw,q.zw);
    float2 d=xx-6*xy+8*yy;
    float2 n=.6*cm(alpha,8*xy-3*xx)+cm(c,d);
    float4 r=float4(n,d); float den=dot(r,r);
    ly=log(max(1e-24,.6*length(alpha)*length(10*xx-48*xy+64*yy)/max(den,1e-30)));
    return r*rsqrt(max(den,1e-30));
}
float3 ramp(float t) { return .51+.46*cos(6.283185307*(float3(0,.105,.19)+t)); }
float3 stretchColor(float ly) { return ramp(option.z+.48+option.w*.20*(ly*.65)/(1+abs(ly*.65)))*(.32+.5*(1-exp(-abs(ly)*1.8))); }
float3 trapColor(float d) { float t=-log(max(1e-12,d));return ramp(option.z+.12+t*.095*option.w)*(.36+.48*(1-exp(-t*.22))); }
float3 inspect(float4 q,float2 c,float2 alpha) {
    float sum=0,count=0,trap=1,ly=0;
    float4 poleA=pt(float2(2,0)),poleB=pt(float2(4,0));
    [loop] for(uint i=0;i<mode.w;i++) {
        q=advance(q,c,alpha,ly);
        if(float(i)>float(mode.w)*.72){sum+=ly;count++;}
        trap=min(trap,min(cd(q,poleA),cd(q,poleB)));
    }
    float period=0;float4 base=q;
    [loop] for(uint j=1;j<=16;j++){q=advance(q,c,alpha,ly);if(period==0&&cd(q,base)<.00002)period=float(j);}
    return float3(sum/max(count,1),period,trap);
}
float3 parameterColor(float2 p) {
    float2 c=mode.x==1?p:param.xy,alpha=mode.x==2?p:param.zw;
    float3 a=0,b=0;
    if(mode.y==3)a=b=inspect(pt(option.xy),c,alpha);
    else {
        if(mode.y!=2)a=inspect(pt(float2(2.4,.8)),c,alpha);
        if(mode.y!=1)b=inspect(pt(float2(2.4,-.8)),c,alpha);
        if(mode.y==1)b=a;if(mode.y==2)a=b;
    }
    float ly=max(a.x,b.x);
    if(mode.z==1)return stretchColor(ly);
    if(mode.z==2)return trapColor(min(a.z,b.z));
    float period=max(a.y,b.y);
    float t=option.z+.08+option.w*(period*.069+log(1+abs(ly))*.10);
    float light=.32+.47*(1-exp(-abs(ly)*1.3));
    if(period<.5){t=option.z+.54+.10*sin(ly*4*option.w);light=.43+.40*clamp(ly,0,1);}
    return ramp(t)*light+float3(.12,.17,.17)*exp(-abs(ly)*14);
}
float3 juliaColor(float2 p) {
    float4 q=pt(p);float sum=0,count=0,ly=0,hit=0,at=0,distance=1,trap=1;
    float4 poleA=pt(float2(2,0)),poleB=pt(float2(4,0));
    [loop] for(uint i=0;i<mode.w;i++) {
        q=advance(q,param.xy,param.zw,ly);if(i>8){sum+=ly;count++;}
        if(mode.z==2)trap=min(trap,min(cd(q,poleA),cd(q,poleB)));
        if(mode.z==0) {
            float da=flags.x!=0?cd(q,refA):2,db=flags.y!=0?cd(q,refB):2;
            float d=min(da,db);
            if(d<.00004){hit=da<=db?1:2;at=float(i);distance=max(d,1e-12);break;}
        }
    }
    if(mode.z==2)return trapColor(trap);
    if(hit>0) {
        float smooth=at+clamp(log(distance/.00004)/min(sum/max(count,1),-.08),-3,3);
        return ramp(option.z+.085+log(1+max(smooth,0))*.105*option.w+(hit-1)*.40)*(.4+.33*(1-exp(-max(smooth,0)*.055)));
    }
    return stretchColor(sum/max(count,1));
}
[numthreads(16,8,1)]
void main(uint3 id:SV_DispatchThreadID) {
    uint x=id.x,y=id.y+sizeRow.z;
    if(x>=sizeRow.x||y>=sizeRow.y||id.y>=sizeRow.w)return;
    float3 col=0;
    [loop] for(uint sy=0;sy<flags.z;sy++) [loop] for(uint sx=0;sx<flags.z;sx++) {
        float2 uv=(float2(x,y)+(float2(sx,sy)+.5)/float(flags.z))/float2(sizeRow.xy);
        float2 p=view.xy+(uv-float2(.5,.5))*float2(view.z,-view.z*view.w);
        col+=mode.x==0?juliaColor(p):parameterColor(p);
    }
    col/=float(flags.z*flags.z);
    uint3 rgb=(uint3)clamp(pow(max(col,0),.88)*255+.5,0,255);
    pixels[uint2(x,y)]=rgb.b|(rgb.g<<8)|(rgb.r<<16);
}
)AMPHLSL";
