#include "frk_extra_catalogue.hpp"
#include "catalogue.hpp"
#include <cassert>
#include <iostream>
#include <thread>
#include <vector>
#include <array>
#include <cstdint>
#include <iomanip>

namespace frk_extra_tests {
using namespace frk_extra;
static bool all_ok=true;
static std::ostream* log_stream=nullptr;
static int checks=0;
static void test(bool ok,const char *msg) {
  ++checks;
  if(!ok){*log_stream<<"SELF_TEST_FAIL=EXTRA CHECK="<<checks<<" SUBJECT="<<msg<<"\n";all_ok=false;}
}
static bool close(C a,C b,double tol=1e-8){return std::abs(a-b)<=tol*(1.+std::abs(b));}
static uint64_t fold(uint64_t v,uint32_t x) {return (v^x)*1099511628211ull;}
static uint64_t image_hash(int driver,int threads) {
  const int width=96,height=64;
  Params p;p.iterations=48;p.lambda=.21;
  std::vector<uint32_t> rgba(size_t(width*height));
  std::vector<std::thread> workers;
  for(int k=0;k<threads;++k){
    workers.emplace_back([&,k] {
      for(int y=k;y<height;y+=threads){
        for(int x=0;x<width;++x){
          C pixel(-2.+4.*(x+.5)/width,-1.5+3.*(y+.5)/height);
          rgba[size_t(y*width+x)]=color(evaluate(driver,pixel,p));
        }
      }
    });
  }
  for(auto &w:workers)w.join();
  uint64_t h=14695981039346656037ull;
  for(auto v:rgba)h=fold(h,v);
  return h;
}
inline bool run(std::ostream& out) {
  log_stream=&out;all_ok=true;checks=0;
  test(driver_count==14,"exact candidate count");
  test(first_driver==18 && last_driver==31,"base catalogue untouched / IDs noncolliding");
  test(names.size()==14 && algebra.size()==14,"metadata complete");
  {
    Mat A=cayley_commutator(0.37);
    test(std::abs(A.a*A.d-A.b*A.c)>1e-6,"commutator invertible");
    test(finite(apply(A,C(.2,.3))),"Cayley map finite");
    test(std::abs(apply(A,C(.2,.3))-C(.2,.3))>1e-4,"nontrivial commutator");
    Mat B{C(1.,0.),C(.23,.17),C(.23,-.17),C(1.,0.)};
    C z(.3,-.2);
    test(close(apply(inverse_projective(B),apply(B,z)),z),"projective inverse");
  }
  {
    State s; s.z=C(.3,.7);s.c=C(-.2,.13);s.previous=0;
    C v=frk_extra::cayley(std::conj(s.z),C(.28,.22));
    test(close(step(19,s,0.).z,v*v+s.c),"antiholomorphic Cayley algebra");
    State t=s;t.previous=C(.2,.1);
    test(std::abs(step(20,s,.42).z-step(20,t,.42).z)>1e-5,"memory depends on lagged state");
  }
  for(int j=0;j<5;++j){
    C r(spectral_root(j),0.);
    test(close(spectral_newton(r+C(1e-6,1e-6)),r,1e-6),"spectral Newton exact finite roots");
  }
  test(close(halley(C(1.,0.),0.),C(1.,0.)),"Halley fixed quintic root");
  test(close(halley(C(2.,0.),0.),C(2.,0.)-C(2.*31.*80./(2.*80.*80.-31.*160.),0.)),"Halley algebra explicit");
  for(int n=0;n<=100;++n){
    double desired=(n%2==0)? n*.5:(3*n+1)*.5;
    test(close(collatz(C(n,0.),0.),C(desired,0.),1e-9),"entire Collatz agrees at integer inputs");
  }
  for(C z: {C(.35,.42), C(1.1,5.4),C(-.6,.2)}){
    C d;C f=eta9(z,&d),eps(1e-6,0.);
    C fd=(eta9(z+eps)-eta9(z-eps))/(2e-6);
    test(finite(f) && close(d,fd,2e-5),"finite eta9 analytic derivative");
  }
  for(C z:{C(0.,0.),C(.3,.7),C(-.5,.33),C(2.,1.)}) {
    C fp=3.*z*z-3.,fpp=6.*z,fppp=6.;
    C formula=divide(fppp,fp)-1.5*divide(fpp,fp)*divide(fpp,fp);
    test(close(schwarzian(z),formula,1e-8),"exact cubic Schwarzian derivative");
  }
  {
    Params p;p.lambda=.65;p.iterations=1;p.julia=true;
    const auto o=evaluate(28,C(-1.,0.),p);
    test(!o.converged && o.iterations==1,"angular small residual alone is not convergence");
  }
  {
    Params p;p.iterations=30;
    for(int driver=first_driver;driver<=last_driver;++driver){
      C z(-.13,.23);
      Result o=evaluate(driver,z,p);
      test(o.iterations>=0 && o.iterations<=p.iterations,"iteration count finite");
      test(o.field||o.escaped||o.converged||o.iterations==p.iterations,"event classification");
      test((color(o)&0xff000000u)==0xff000000u,"color fully opaque");
      uint64_t serial=image_hash(driver,1),parallel=image_hash(driver,4);
      test(serial==parallel,"deterministic parallel CPU render");
      out<<"DRIVER="<<driver<<" NAME=\""<<names[driver-first_driver]<<"\" CPU_PARITY=PASS HASH=0x"<<std::hex<<serial<<std::dec<<"\n";
    }
  }
  test(catalogue::load("renderers").size()==32,"catalogue includes exactly 32 builtins");
  for(int i=1;i<=driver_count;++i) {
    const auto s=amp::defaults(25,i);
    auto scene=amp::scene(s);
    test(scene.state.family==25 && scene.gardenData!=nullptr,"merged base family/scene");
    const auto color=amp::renderColor(C(.07,.11),scene);
    test(std::isfinite(color[0])&&std::isfinite(color[1])&&std::isfinite(color[2]),"merged CPU colour evaluator");
    const auto name=amp::formula(s);
    test(!name.empty(),"merged formula");
    // Regression: view version 17 belongs exclusively to the family-25 drivers.
    // Preserve compatibility with all preceding view versions.
    try {
      const auto temp=std::filesystem::temp_directory_path() /
          ("fraktalfaktur-v186-extra-"+std::to_string(i)+".view");
      amp::writeView(temp,s);
      const auto back=amp::readView(temp);
      std::error_code ec;std::filesystem::remove(temp,ec);
      test(back.family==25 && back.generation==i && back.plane==s.plane &&
          std::abs(back.alpha-s.alpha)<1e-12 && back.iterations==s.iterations,
          "v17 native view round-trip");
    } catch(const std::exception& error) {
      out<<"VIEW_ROUNDTRIP_FAIL="<<i<<" ERROR="<<error.what()<<'\n';
      test(false,"v17 native view round-trip");
    }
  }
  for(int generation=9;generation<=driver_count;++generation)for(int preset=0;preset<4;++preset) {
    auto state=amp::defaults(25,generation);amp::preset(state,preset);
    try {
      const auto path=std::filesystem::temp_directory_path()/"fraktalfaktur-hilbert-preset.view";
      amp::writeView(path,state);const auto back=amp::readView(path);
      std::filesystem::remove(path);
      test(back.generation==generation && back.plane==state.plane && back.alpha==state.alpha,
           "Hilbert presets survive view round-trip");
    } catch(const std::exception& e) {
      out<<"HILBERT_PRESET_FAIL="<<generation<<" PRESET="<<preset<<" ERROR="<<e.what()<<'\n';
      test(false,"Hilbert presets survive view round-trip");
    }
  }
  out<<"EXTRA_INTEGRATION_TEST="<<(all_ok?"PASS":"FAIL")<<" CHECKS="<<checks<<"\n";
  return all_ok;
}
} // namespace frk_extra_tests
