#pragma once
#include "frk_extra_math.hpp"
// Integration metadata for existing v1.8 catalogue. The host remains the source
// of truth for UI, input bindings, shader dispatch, CPU/GPU work scheduling.
namespace frk_extra {
struct Entry {
  int driver_id;
  const char *label;
  double center_real, center_imag, vertical_span;
  double default_lambda;
  int suggested_iterations;
  bool start_at_pixel;
  bool supports_julia_plane;
};
static constexpr std::array<Entry,driver_count> extra_catalogue={{
  {18,"Cayley commutator", 2.2, 1.1, 3.3,.30,160,false,true},
  {19,"Tricorn Cayley",2.1,-1.4,3.6,.30,160,false,true},
  {20,"GP echo memory",0.,0.,4.7,.50,160,false,true},
  {21,"Spectral Newton web",-3.15,0.,5.5,.30,90,true,false},
  {22,"Halley quintic web",0.,0.,3.9,.45,90,true,false},
  {23,"Entire Collatz parity",0.,0.,7.5,.0,100,true,false},
  {24,"Dirichlet eta Newton",.15,11.5,20.,.0,90,true,false},
  {25,"Schwarzian pole web",0.,0.,5.6,.13,160,false,true},
  {26,"QQRH Hilbert overlap",0.,0.,3.,.65,96,true,false},
  {27,"QQRH Hardy shift fractal",0.,0.,3.6,.65,160,false,true},
  {28,"QQRH angular kernel Newton",0.,0.,4.,.65,96,true,false},
  {29,"QQRH projective metric",0.,0.,3.,.65,96,true,false},
  {30,"QQRH overlap kernel Newton",-.6,0.,5.,.65,120,true,false},
  {31,"QQRH kernel feedback fractal",0.,0.,3.6,.65,160,false,true}
}};
}
