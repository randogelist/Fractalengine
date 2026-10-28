#pragma once
static const char gardenShader[] = R"GARDENHLSL(
cbuffer Settings : register(b0) {
    uint4 sizeRow;
    uint4 mode;
    float4 view;
    float4 param;
    float4 option;
    float4 refA;
    float4 refB;
    uint4 flags;
    uint4 garden;    // family, generation, method, degree
    uint4 record;    // record count
};
Texture2D<float4> coefficients : register(t0);
RWTexture2D<uint> pixels : register(u0);
float2 cm(float2 a,float2 b){return float2(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);}
float2 divide(float2 a,float2 b){float scale=max(abs(b.x),abs(b.y));float2 q=b/scale;return float2(dot(a,q),a.y*q.x-a.x*q.y)/(scale*dot(q,q));}
float2 power(float2 z,uint n){float2 p=float2(1,0);[loop]while(n>0){if((n&1)!=0)p=cm(p,z);z=cm(z,z);n>>=1;}return p;}
struct Evaluation {float2 f;float2 first;float2 second;float2 base;float magnitude;};
struct Jet {float2 f;float2 first;float2 second;};
Jet jetAdd(Jet a,Jet b){Jet v;v.f=a.f+b.f;v.first=a.first+b.first;v.second=a.second+b.second;return v;}
Jet jetMul(Jet a,Jet b){Jet v;v.f=cm(a.f,b.f);v.first=cm(a.first,b.f)+cm(a.f,b.first);v.second=cm(a.second,b.f)+2*cm(a.first,b.first)+cm(a.f,b.second);return v;}
Jet jetScale(Jet a,float s){a.f*=s;a.first*=s;a.second*=s;return a;}
Evaluation evaluate(float2 z,float2 a,float2 b){
    Evaluation e;e.f=0;e.first=0;e.second=0;e.base=float2(1,0);e.magnitude=0;
    if(garden.x==5){
        Jet p,q;p.f=z;p.first=float2(1,0);p.second=0;q.f=float2(1,0);q.first=0;q.second=0;
        [loop]for(uint k=0;k<garden.y;k++){
            Jet sum=jetAdd(p,q),nextP=jetScale(jetMul(p,q),4),nextQ=jetMul(sum,sum);
            float norm=max(max(length(nextP.f),length(nextP.first)),max(length(nextP.second),max(length(nextQ.f),max(length(nextQ.first),length(nextQ.second)))));
            if(!(norm>0)||!(abs(norm)<=3.0e38))return e;
            // One common scale preserves the local polynomial derivative jets.
            p=jetScale(nextP,1/norm);q=jetScale(nextQ,1/norm);
        }
        e.f=p.f-cm(a,q.f);e.first=p.first-cm(a,q.first);e.second=p.second-cm(a,q.second);
        e.magnitude=length(p.f)+(1+length(a))*length(q.f);
    }else if(garden.x==1){
        uint degree=garden.w;if(dot(b,b)==0)degree=dot(a,a)==0?1:(1u<<garden.y);
        bool reverse=dot(z,z)>1;float2 w=reverse?divide(float2(1,0),z):z;
        bool tiny=dot(z,z)<1e-12;
        float2 ap[5],bp[5];ap[0]=bp[0]=float2(1,0);
        [loop]for(uint k=1;k<=garden.y;k++){ap[k]=cm(ap[k-1],a);bp[k]=cm(bp[k-1],b);}
        [loop]for(uint i=0;i<record.x;i++){
            float4 r=coefficients.Load(int3(i,0,0));uint exponent=(uint)r.x;
            float2 coefficient=r.w*cm(ap[(uint)r.y],bp[(uint)r.z]);if(dot(coefficient,coefficient)==0)continue;
            uint k=reverse?degree-exponent:exponent;float2 value=cm(coefficient,power(w,k));e.f+=value;e.magnitude+=length(value);
            if(tiny){if(exponent>0)e.first+=float(exponent)*cm(coefficient,power(z,exponent-1));if(exponent>1)e.second+=float(exponent*(exponent-1))*cm(coefficient,power(z,exponent-2));}
            else{e.first+=float(exponent)*value;e.second+=float(exponent*(exponent-1))*value;}
        }
        e.base=tiny?float2(1,0):z;
    }else if(garden.x==4){
        bool reverse=dot(z,z)>1,tiny=dot(z,z)<1e-12;float2 w=reverse?divide(float2(1,0),z):z;
        [loop]for(int i=(int)garden.w;i>=0;i--){
            int k=reverse?(int)garden.w-i:i;float4 r=coefficients.Load(int3(k,0,0));float2 coefficient=float2(r.x,0)-a*r.y;
            e.f=cm(e.f,w)+coefficient;e.magnitude=e.magnitude*length(w)+length(coefficient);
            if(tiny){if(k>0)e.first+=float(k)*cm(coefficient,power(z,k-1));if(k>1)e.second+=float(k*(k-1))*cm(coefficient,power(z,k-2));}
            else{e.first=cm(e.first,w)+float(k)*coefficient;e.second=cm(e.second,w)+float(k*(k-1))*coefficient;}
        }
        e.base=tiny?float2(1,0):z;
    }else{
        int m=(int)garden.w;float2 oneMinus=float2(1,0)-z;bool reverse=dot(z,z)>dot(oneMinus,oneMinus);
        e.base=reverse?z:oneMinus;float2 w=reverse?divide(oneMinus,z):divide(z,oneMinus);
        [loop]for(int i=m;i>=0;i--){int k=reverse?m-i:i;float4 r=coefficients.Load(int3(k,0,0));float2 value=float2(r.x,0)-a*r.y;e.f=cm(e.f,w)+value;e.magnitude=e.magnitude*length(w)+length(value);}
        [loop]for(int i=m-1;i>=0;i--){int k=reverse?m-1-i:i;e.first=cm(e.first,w)+float2(coefficients.Load(int3(k,0,0)).z,0);}
        [loop]for(int i=m-2;i>=0;i--){int k=reverse?m-2-i:i;e.second=cm(e.second,w)+float2(coefficients.Load(int3(k,0,0)).w,0);}
    }
    return e;
}
float phase(float2 p){
    float ax=abs(p.x),ay=abs(p.y),t=min(ax,ay)/max(max(ax,ay),1e-30),q=t*t;
    float a=t*(.9998660+q*(-.3302995+q*(.1801410+q*(-.0851330+q*.0208351))));
    if(ay>ax)a=1.57079632679-a;if(p.x<0)a=3.14159265359-a;if(p.y<0)a=-a;return a;
}
float3 rainbow(float t){return .52+.46*cos(6.283185307*(float3(0,.33,.67)+t));}
float3 crownColour(float2 p){
    float2 z=mode.x==0?p:option.xy,q=mode.x==1?p:param.xy,b=mode.x==2?p:param.zw;
    float2 roots[81];
    [loop]for(uint k=0;k<record.x;k++){
        float4 r=coefficients.Load(int3(k,0,0));
        float2 root=r.xy;
        if(mode.x==1){
            root=0;uint residue=(uint)r.z;
            [loop]for(int j=(int)r.w;j>=1;j--)root=cm(root,q)+float2(2*float((residue>>j)&1u)-1,0);
        }
        if(!all(abs(root)<=3.0e38))return float3(.09,.10,.13);
        roots[k]=root;
    }
    uint status=1,steps=0;float error=1,label=0;
    [loop]for(uint i=0;i<mode.w;i++){
        if(!all(abs(z)<=3.0e38)||length(z)>1e7){status=3;break;}
        steps=i;error=3.0e38;
        [loop]for(uint k=0;k<record.x;k++){
            float distance=length(z-roots[k])/(1+length(roots[k]));
            if(distance<error){error=distance;label=(float(k)+.5)/float(record.x);}
        }
        if(error<3e-6){status=0;break;}
        float2 first=0,second=0;
        [loop]for(uint k=0;k<record.x;k++){
            float2 inv=divide(float2(1,0),z-roots[k]);first+=inv;second+=cm(inv,inv);
        }
        float2 denominator=garden.z==1?cm(first,first)+second:first;
        if(dot(denominator,denominator)<1e-26){status=2;break;}
        float2 delta=cm(b,divide(garden.z==1?2*first:float2(1,0),denominator));
        if(!all(abs(delta)<=3.0e38)){status=2;break;}
        z-=delta;steps=i+1;
    }
    if(status==2)return float3(.19,.095,.065);
    if(status==3)return float3(.09,.10,.13);
    if(status!=0)return rainbow(option.z+.72+.015*log(1+error))*.12;
    float t=option.z+label,light=.37+.44*exp(-float(steps)*.035);
    if(mode.z==0)light*=.84+.16*cos(log(1+float(steps))*option.w*3.8);
    if(mode.z==1){t=option.z+.04+option.w*log(1+float(steps))*.13;light=.78;}
    if(mode.z==2){t=option.z+.08-option.w*(log(max(error,3e-7))/2.302585093)*.065;light=.72;}
    return rainbow(t)*light;
}
float3 gardenColour(float2 p){
    if(garden.x==3)return crownColour(p);
    float2 z=mode.x==0?p:option.xy,a=mode.x==1?p:param.xy,b=mode.x==2?p:param.zw;
    uint status=1,steps=0;float error=1;
    [loop]for(uint i=0;i<mode.w;i++){
        if(!all(abs(z)<=3.0e38)||length(z)>1e7){status=3;break;}
        if(garden.x==4&&a.x==0&&a.y==0){
            steps=i;error=length(z)/(1+length(z));
            if(error<3e-6){z=0;status=0;break;}
            float2 delta=cm(b,z)*(garden.z==1?2./float(garden.w+1):1./float(garden.w));
            if(!all(abs(delta)<=3.0e38)){status=2;break;}
            z-=delta;steps=i+1;continue;
        }
        Evaluation e=evaluate(z,a,b);steps=i;
        if(!(e.magnitude>0)||!(abs(e.magnitude)<=3.0e38)){status=2;break;}
        error=length(e.f)/e.magnitude;
        if(garden.x==5&&all(e.f==float2(0,0))){status=0;break;}
        float norm=max(length(e.f),max(length(e.first),length(e.second)));
        if(!(norm>0)||!(abs(norm)<=3.0e38)){status=2;break;}
        e.f/=norm;e.first/=norm;e.second/=norm;
        float2 denominator=garden.z==1?2*cm(e.first,e.first)-cm(e.f,e.second):e.first;
        if(dot(denominator,denominator)<1e-26){status=error<3e-6?0:2;break;}
        float2 delta=cm(e.base,divide(garden.z==1?2*cm(e.f,e.first):e.f,denominator));
        if(garden.x==2||garden.x==4||garden.x==5)delta=cm(delta,b);
        if(!all(abs(delta)<=3.0e38)){status=2;break;}
        if(error<3e-6){status=0;break;}
        z-=delta;steps=i+1;
    }
    if(status==2)return float3(.19,.095,.065);
    if(status==3)return float3(.09,.10,.13);
    if(status!=0)return rainbow(option.z+.72+.015*log(1+error))*.12;
    float2 root=z-float2(garden.x==2?.5:0,0);
    float angle=phase(root)/6.283185307;
    float t=option.z+angle+.06*log(1+length(root));float light=.37+.44*exp(-float(steps)*.035);
    if(mode.z==0)light*=.84+.16*cos(log(1+float(steps))*option.w*3.8);
    if(mode.z==1){t=option.z+.04+option.w*log(1+float(steps))*.13;light=.78;}
    if(mode.z==2){t=option.z+.08-option.w*(log(max(error,garden.x>=3?3e-7:1e-16))/2.302585093)*.065;light=.72;}
    return rainbow(t)*light;
}
[numthreads(16,8,1)]
void main(uint3 id:SV_DispatchThreadID){
    uint x=id.x,y=id.y+sizeRow.z;if(x>=sizeRow.x||y>=sizeRow.y||id.y>=sizeRow.w)return;
    float3 col=0;
    [loop]for(uint sy=0;sy<flags.z;sy++)[loop]for(uint sx=0;sx<flags.z;sx++){
        float2 uv=(float2(x,y)+(float2(sx,sy)+.5)/float(flags.z))/float2(sizeRow.xy);
        float2 p=view.xy+(uv-.5)*float2(view.z,-view.z*view.w);col+=gardenColour(p);
    }
    col/=float(flags.z*flags.z);uint3 rgb=(uint3)clamp(pow(max(col,0),.88)*255+.5,0,255);
    pixels[uint2(x,y)]=rgb.b|(rgb.g<<8)|(rgb.r<<16);
}
)GARDENHLSL";
