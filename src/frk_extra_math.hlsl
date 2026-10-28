// FRAKTALFAKTUR v1.8 additional candidate GPU evaluation functions.
// Shader Model 5.0 / D3D11. INCLUDE-ONLY; does not take over the base renderer.
// The host's pixel/compute shader must call frkx_evaluate(driver, pixel, lambda, julia, juliaC)
// for driver IDs 18..31 and use frkx_color() or the host's own existing palette.
#ifndef FRAKTALFAKTUR_EXTRA_CANDIDATES_V181
#define FRAKTALFAKTUR_EXTRA_CANDIDATES_V181
static const float FRKX_PI = 3.14159265358979323846;
float2 frkx_cmul(float2 a,float2 b) {return float2(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);}
float2 frkx_cconj(float2 a) {return float2(a.x,-a.y);}
float2 frkx_cdiv(float2 a,float2 b) {
  float d=dot(b,b);
  if(!(d>1e-25) || !(d<1e26)) return float2(1e12,1e12);
  return float2(dot(a,b),a.y*b.x-a.x*b.y)/d;
}
float2 frkx_cexp(float2 z) {
  float a=exp(clamp(z.x,-50.0,50.0));
  return a*float2(cos(z.y),sin(z.y));
}
float2 frkx_csin(float2 z){
  float hi=exp(clamp(z.y,-15.0,15.0));
  float lo=exp(clamp(-z.y,-15.0,15.0));
  return .5*float2(sin(z.x)*(hi+lo),cos(z.x)*(hi-lo));
}
float2 frkx_apply(float2 A,float2 B,float2 C,float2 D,float2 z) {
  return frkx_cdiv(frkx_cmul(A,z)+B,frkx_cmul(C,z)+D);
}
struct FrkxMat {float2 a,b,c,d;};
FrkxMat frkx_matmul(FrkxMat X,FrkxMat Y) {
  FrkxMat R;
  R.a=frkx_cmul(X.a,Y.a)+frkx_cmul(X.b,Y.c);
  R.b=frkx_cmul(X.a,Y.b)+frkx_cmul(X.b,Y.d);
  R.c=frkx_cmul(X.c,Y.a)+frkx_cmul(X.d,Y.c);
  R.d=frkx_cmul(X.c,Y.b)+frkx_cmul(X.d,Y.d);
  return R;
}
FrkxMat frkx_matinverse(FrkxMat A) {
  FrkxMat R;
  R.a=A.d;R.b=-A.b;R.c=-A.c;R.d=A.a;
  return R;
}
FrkxMat frkx_commutator(float lambda) {
  float2 a=float2(.48,.13)*(1.+.20*lambda);
  float2 b=float2(-.26,.43)*(1.+.15*lambda);
  FrkxMat A,B;
  A.a=float2(1.,0.);A.b=a;A.c=frkx_cconj(a);A.d=float2(1.,0.);
  B.a=float2(1.,0.);B.b=b;B.c=frkx_cconj(b);B.d=float2(1.,0.);
  FrkxMat R=frkx_matmul(frkx_matmul(frkx_matmul(A,B),frkx_matinverse(A)),frkx_matinverse(B));
  float scale=max(1e-12,max(max(length(R.a),length(R.b)),max(length(R.c),length(R.d))));
  R.a/=scale;R.b/=scale;R.c/=scale;R.d/=scale;
  return R;
}
float2 frkx_cayley(float2 z,float2 a){return frkx_cdiv(z-a,float2(1.,0.)-frkx_cmul(frkx_cconj(a),z));}
float frkx_spectral_root(int i) {
  if(i==0)return -1.;
  if(i==1)return -2.21194826822;
  if(i==2)return -3.13099540926;
  if(i==3)return -4.63321978045;
  return -5.42927980568;
}
void frkx_eta9(float2 z,out float2 f,out float2 d){
  f=float2(0.,0.);d=float2(0.,0.);
  float sign=1.;
  for(int j=1;j<=9;++j){
    float ln=log((float)j);
    float2 t=sign*frkx_cexp(-z*ln);
    f+=t;d-=ln*t;sign=-sign;
  }
}
float2 frkx_schwarzian(float2 z) {
  float2 z2=frkx_cmul(z,z),t=z2-float2(1.,0.);
  return frkx_cdiv(-4.*z2-float2(2.,0.),frkx_cmul(t,t));
}
// Finite Hardy model: T(z)=Theta(200*z)/Theta(0), degree 12.
// These are the unchanged floating-point coefficients in theta_coefficients.hpp.
// No infinite-dimensional positivity or RH certification is implied.
static const float FRKH_B[13]={
  1.00000000000000000,4.62099862308439224,9.93336215156702629,
  13.3948210238812724,12.8491158736715363,9.40843442569675048,
  5.50310055257359476,2.65455036910845843,1.08130879752268716,
  .378823079526009854,.115823084088392592,.0312759502906683751,
  .00753360834543220954
};
void frkh_model(float radius,out float v[13],out float correlation[13],out float angularScale){
  float R=clamp(radius,.1,1.25),power=1.,norm2=0.;
  [unroll]for(int i=0;i<13;++i){v[i]=FRKH_B[i]*power;norm2+=v[i]*v[i];power*=R;}
  float invNorm=rsqrt(norm2);
  [unroll]for(int j=0;j<13;++j)v[j]*=invNorm;
  angularScale=0.;
  [unroll]for(int k=0;k<13;++k){
    float sum=0.;
    [unroll]for(int n=0;n<13;++n){if(n+k<13)sum+=v[n]*v[n+k];}
    correlation[k]=sum;
    angularScale+=float(k*k)*sum;
  }
}
void frkh_correlation(float2 z,float correlation[13],out float2 f,out float2 d){
  // Simultaneous Horner value/derivative, including the zero constant term.
  f=float2(correlation[12]/correlation[1],0.);d=float2(0.,0.);
  [unroll]for(int k=11;k>=0;--k){
    d=frkx_cmul(d,z)+f;
    f=frkx_cmul(f,z)+float2(k==0 ? 0. : correlation[k]/correlation[1],0.);
  }
}
void frkh_angular(float2 z,float correlation[13],float angularScale,out float2 f,out float2 d){
  float2 previous=float2(0.,0.),current=float2(1.,0.);
  float2 previousD=float2(0.,0.),currentD=float2(0.,0.);
  f=float2(0.,0.);d=float2(0.,0.);
  [unroll]for(int k=1;k<=12;++k){
    float weight=float(k)*correlation[k]/angularScale;
    f+=weight*current;d+=weight*currentD;
    float2 next=2.*frkx_cmul(z,current)-previous;
    float2 nextD=2.*current+2.*frkx_cmul(z,currentD)-previousD;
    previous=current;current=next;previousD=currentD;currentD=nextD;
  }
}
float frkh_coherence(float2 z,float v[13]){
  // Rescale all feature coordinates by z^-12 outside the unit disk.
  // The common scalar cancels, keeping distant zooms free of power overflow.
  float largest=max(abs(z.x),abs(z.y));
  if(!(largest<1e30))return 0.;
  bool reverse=largest>1. || dot(z,z)>1.;
  float2 base=z;
  if(reverse){
    float2 scaled=z/max(largest,1e-20);
    base=float2(scaled.x,-scaled.y)/(max(largest,1e-20)*dot(scaled,scaled));
  }
  float2 power=float2(1.,0.),overlap=float2(0.,0.);
  float norm2=0.,anchorNorm2=0.;
  [unroll]for(int i=0;i<13;++i){
    int n=reverse ? 12-i : i;
    float2 feature=v[n]*power;
    overlap+=v[n]*feature;norm2+=dot(feature,feature);anchorNorm2+=v[n]*v[n];
    power=frkx_cmul(power,base);
  }
  return norm2>0. ? saturate(dot(overlap,overlap)/(norm2*anchorNorm2)) : 0.;
}
// H(z)=<w(z),w(1)>=sum v_n^2 z^n, with an exact polynomial jet.
void frkh_kernel(float2 z,float v[13],out float2 f,out float2 d){
  f=float2(v[12]*v[12],0.);d=float2(0.,0.);
  [unroll]for(int k=11;k>=0;--k){
    d=frkx_cmul(d,z)+f;
    f=frkx_cmul(f,z)+float2(v[k]*v[k],0.);
  }
}
float frkh_degree_variance(float t,float v[13],bool reverse){
  // Positive pair formula for the degree variance. Feature scaling prevents
  // the squared tiny leading coefficient from underflowing in reversed charts.
  int first=reverse ? 12 : 0;
  int second=reverse ? 11 : 1;
  if(t<1e-12){float ratio=v[second]/v[first];return (ratio*ratio)*t;}
  float feature[13],weights[13];
  float power=1.,largest=0.,base=sqrt(t);
  [unroll]for(int i=0;i<13;++i){
    int index=reverse ? 12-i : i;
    feature[i]=v[index]*power;
    largest=max(largest,abs(feature[i]));power*=base;
  }
  float total=0.;
  [unroll]for(int j=0;j<13;++j){
    float scaled=feature[j]/max(largest,1e-30);
    weights[j]=scaled*scaled;total+=weights[j];
  }
  float pairs=0.;
  [unroll]for(int a=0;a<13;++a){
    [unroll]for(int b=0;b<13;++b){
      if(b>a)pairs+=float((b-a)*(b-a))*weights[a]*weights[b];
    }
  }
  return pairs/(total*total);
}
float frkh_projective_metric(float2 z,float v[13]){
  float largest=max(abs(z.x),abs(z.y));
  if(!(largest<1e30))return 0.;
  if(largest<=1.){
    float t=dot(z,z);
    if(t<1e-12){float ratio=v[1]/v[0];return ratio*ratio;}
    if(t<=1.)return frkh_degree_variance(t,v,false)/t;
  }
  float2 scaled=z/largest;
  float inverse=1./largest;
  float u=(inverse*inverse)/dot(scaled,scaled);
  return u*frkh_degree_variance(u,v,true);
}
struct FrkxResult {
  int iterations,basin,escaped,converged;
  float2 z;
  float trap;
  int field;
  float value;
};
bool frkx_seed_driver(int driver) {return (driver>=21 && driver<=24) || driver==26 || driver==28 || driver==29 || driver==30;}
bool frkx_root_driver(int driver) {return driver==21 || driver==22 || driver==24 || driver==28 || driver==30;}
FrkxResult frkx_evaluate(int driver,float2 pixel,float lambda,bool julia,float2 juliaC,int maxIter){
  FrkxResult r;
  r.iterations=0;r.basin=-1;r.escaped=0;r.converged=0;r.trap=1e12;r.field=0;r.value=0.;
  float hilbertV[13],hilbertC[13];float angularScale=1.;
  [branch]if(driver>=26){
    frkh_model(lambda,hilbertV,hilbertC,angularScale);
    if(driver==26){r.field=1;r.value=frkh_coherence(pixel,hilbertV);r.z=pixel;return r;}
    if(driver==29){
      float metric=frkh_projective_metric(pixel,hilbertV);
      r.field=1;r.value=metric/(1.+metric);r.z=pixel;return r;
    }
  }
  bool seed=frkx_seed_driver(driver);
  float2 c=seed ? float2(0.,0.) : (julia ? juliaC : pixel);
  float2 z=(seed || julia) ? pixel : float2(0.,0.);
  float2 prev=float2(0.,0.);
  float bailout=(driver==23)?72.:(frkx_root_driver(driver)?20000.:64.);
  FrkxMat comm=frkx_commutator(lambda);
  [loop] for(int n=0;n<1024;++n){
    if(n>=maxIter)break;
    if(!(dot(z,z) < bailout*bailout) || !all(abs(z)<1e10) ||
        (driver==23 && abs(z.y)>16.) ||
        (driver==24 && (z.x < -24. || z.x > 28.))){r.escaped=1;r.iterations=n;break;}
    r.trap=min(r.trap,length(z));
    float2 hilbertF=float2(0.,0.),hilbertD=float2(0.,0.);
    if(driver==21){
      float best=1e10;int nearest=-1;
      [unroll]for(int k=0;k<5;++k){
        float t=length(z-float2(frkx_spectral_root(k),0.));
        if(t<best){best=t;nearest=k;}
      }
      if(best<1e-4){r.converged=1;r.basin=nearest;r.iterations=n;break;}
    }else if(driver==22){
      float2 z2=frkx_cmul(z,z),z4=frkx_cmul(z2,z2);
      float2 f=frkx_cmul(z4,z)-float2(1.,0.)+(.30*lambda)*z2;
      if(length(f)<1e-4){r.converged=1;r.basin=(int)floor((atan2(z.y,z.x)+FRKX_PI)/(2.*FRKX_PI)*5.);r.iterations=n;break;}
    }else if(driver==24){
      float2 f,d;frkx_eta9(z,f,d);
      if(length(f)<1e-4){r.converged=1;r.basin=(int)floor((atan2(z.y,z.x)+FRKX_PI)/(2.*FRKX_PI)*8.);r.iterations=n;break;}
    }else if(driver==28 || driver==30){
      if(driver==28)frkh_angular(z,hilbertC,angularScale,hilbertF,hilbertD);
      else frkh_kernel(z,hilbertV,hilbertF,hilbertD);
      if(dot(hilbertD,hilbertD)>1e-25 && length(hilbertF)<1e-4 &&
         length(frkx_cdiv(hilbertF,hilbertD))<1e-4*(1.+length(z))){r.converged=1;r.basin=driver==30 ? -2 : (int)floor((atan2(z.y,z.x)+FRKX_PI)/(2.*FRKX_PI)*12.);r.iterations=n;break;}
    }
    float2 nextZ=z;
    if(driver==18){
      float2 w=frkx_apply(comm.a,comm.b,comm.c,comm.d,z);
      nextZ=frkx_cmul(w,w)+c;
    }else if(driver==19){
      float2 a=float2(.28,.22)*(1.+.17*lambda);
      float2 w=frkx_cayley(frkx_cconj(z),a);
      nextZ=frkx_cmul(w,w)+c+(.14*lambda)*w;
    }else if(driver==20){
      nextZ=frkx_cmul(z,z)+c+frkx_cmul(float2(lambda,.22*lambda),z-prev);
      prev=z;
    }else if(driver==21){
      float2 s=float2(0.,0.);
      [unroll]for(int k=0;k<5;++k){s+=frkx_cdiv(float2(1.,0.),z-float2(frkx_spectral_root(k),0.));}
      nextZ=z-frkx_cdiv(float2(1.,0.),s);
    }else if(driver==22){
      float2 z2=frkx_cmul(z,z),z3=frkx_cmul(z2,z),z4=frkx_cmul(z2,z2);
      float2 P=frkx_cmul(z4,z)-float2(1.,0.)+(.30*lambda)*z2;
      float2 d=5.*z4+(.60*lambda)*z;
      float2 dd=20.*z3+float2(.60*lambda,0.);
      nextZ=z-frkx_cdiv(2.*frkx_cmul(P,d),2.*frkx_cmul(d,d)-frkx_cmul(P,dd));
    }else if(driver==23){
      float2 s=frkx_csin(z*(FRKX_PI*.5));
      nextZ=z*.5+frkx_cmul(z+float2(.5,0.),frkx_cmul(s,s))+
           (.125*lambda)*frkx_csin(z*FRKX_PI);
    }else if(driver==24){
      float2 f,d;frkx_eta9(z,f,d);nextZ=z-frkx_cdiv(f,d);
    }else if(driver==25){
      nextZ=frkx_cmul(z,z)+c+lambda*frkx_schwarzian(z);
    }else if(driver==27){
      float2 f,d;frkh_correlation(z,hilbertC,f,d);
      nextZ=frkx_cmul(z,z)+c-.35*frkx_cdiv(f,d);
    }else if(driver==28 || driver==30){
      nextZ=z-frkx_cdiv(hilbertF,hilbertD);
    }else if(driver==31){
      float2 f,d;frkh_kernel(z,hilbertV,f,d);
      nextZ=frkx_cmul(z,z)+c-.35*frkx_cdiv(f,d);
    }
    z=nextZ;
    r.iterations=n+1;
  }
  r.z=z;
  return r;
}
float4 frkx_color(FrkxResult r,int palette){
  if(r.field!=0){
    float q=saturate(r.value);
    float3 rgb=float3(.035+.90*q,.08+.65*sqrt(q),.16+.75*(1.-q));
    rgb*=.78+.22*cos(12.*FRKX_PI*q);
    if(palette==1)rgb=rgb.yzx;
    if(palette==2)rgb=rgb.zyx;
    return float4(saturate(rgb),1.);
  }
  if(r.converged==0 && r.escaped==0)return float4(.0235,.0745,.1294,1.);
  float phase=atan2(r.z.y,r.z.x);
  float magn=length(r.z);
  float t=(r.converged!=0)?(.173*r.basin+.009*r.iterations+.05*phase):
       (.023*r.iterations+.16*phase/(2.*FRKX_PI)+.012*log(1.+magn));
  if(r.converged!=0 && r.basin==-2)
    t=.47*log(1.+magn)+.10*sin(phase)+.006*r.iterations;
  if(palette==1)t+=.19;
  if(palette==2)t=-t;
  float bright=(r.converged!=0)?(.69+.25*exp(-r.iterations/25.)):
      (.42+.52*(1.-exp(-r.iterations/30.)));
  float3 rgb=bright*(.54+.46*cos(2.*FRKX_PI*t+float3(.1,2.1,4.2)));
  return float4(saturate(rgb),1.);
}
#endif
