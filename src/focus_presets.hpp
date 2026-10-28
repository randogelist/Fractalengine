#pragma once
#include "amp_core.hpp"

// Presets change the experimental dynamics while retaining the selected RH driver.
// Keeping them here also lets the native test runner exercise exactly the UI values.
namespace focus_presets {
constexpr int BankCount=3;
constexpr int PresetsPerBank=4;

inline const char *bankName(int bank) {
    static const char *names[]={"Original Lab","Currents","Weaves"};
    return names[std::clamp(bank,0,BankCount-1)];
}

inline const char *presetName(int bank,int index) {
    static const char *names[BankCount][PresetsPerBank]={
        {"Classic","Newton basins","Tricorn weave","Cayley currents"},
        {"Cayley currents","Conjugate currents","Alternating tides","Inertial eddies"},
        {"Tricorn weave","Braided tricorn","Twisted currents","Cayley weave"}
    };
    return names[std::clamp(bank,0,BankCount-1)][std::clamp(index,0,PresetsPerBank-1)];
}

inline void apply(amp::State &state,int bank,int index) {
    bank=std::clamp(bank,0,BankCount-1);
    index=std::clamp(index,0,PresetsPerBank-1);
    state.rhExperiment={};
    state.alpha=amp::C(.45,0);
    state.c=amp::defaults(20,state.generation).c;
    auto &e=state.rhExperiment;
    auto tricorn=[&]() {
        e.map=2;state.c=amp::C(-.12,.72);state.alpha=amp::C(.18,0);
        e.coloring=3;e.frequency=2.5;e.trapRadius=.6;e.trapAngle=.4;
    };
    auto cayley=[&]() {
        e.map=5;state.c=amp::C(-.56,.23);state.alpha=amp::C(.16,0);
        e.coloring=5;e.frequency=3;
    };
    if(bank==0) {
        if(index==1) {e.map=1;state.c=amp::C{};state.alpha=amp::C(1,0);e.coloring=1;}
        else if(index==2)tricorn();
        else if(index==3)cayley();
    } else if(bank==1) {
        if(index==0)cayley();
        else if(index==1) {
            e.map=6;state.c=amp::C(-.12,.72);state.alpha=amp::C(.12,0);
            e.coloring=5;e.frequency=3;
        } else if(index==2) {
            e.map=7;state.c=amp::C(-.4,.23);state.alpha=amp::C(.10,0);
            e.coloring=4;e.frequency=1.7;e.relaxation=.85;
        } else {
            e.map=8;state.c=amp::C(-.56,.23);state.alpha=amp::C(.12,0);
            e.coloring=2;e.trapRadius=.55;e.memory=.24;e.relaxation=.9;
        }
    } else {
        if(index==0)tricorn();
        else if(index==1) {
            e.map=10;e.blend=.78;state.c=amp::C(-.12,.72);state.alpha=amp::C(.16,0);
            e.coloring=3;e.trapAngle=.6;e.frequency=2.5;
        } else if(index==2) {
            cayley();e.map=9;e.twist=.7;state.alpha=amp::C(.14,0);
        } else {
            e.map=11;e.blend=.85;e.twist=-.45;state.c=amp::C(-.12,.72);
            state.alpha=amp::C(.11,0);e.coloring=3;e.frequency=2.8;e.trapAngle=.4;
        }
    }
    state.plane=1;
    state.cameras[1]=amp::Camera{0,0,4};
    state.preset=index;
}
} // namespace focus_presets
