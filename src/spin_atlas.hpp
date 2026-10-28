#pragma once
#include "spin_math.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <vector>

namespace spin {
struct AtlasStats {
    uint64_t cells = 0;
    int gridWidth = 0, gridHeight = 0, burnIn = 0, window = 0, startsPerCell = 1;
    double minLyapunov = 0, maxLyapunov = 0, meanLyapunov = 0, meanAreaRate = 0;
};

// Every parameter cell starts from the SAME seeded initial point and tangent.
// This estimates finite-time tangent growth for one initial condition; it is
// neither an exhaustive basin survey nor a certified asymptotic exponent.
class AtlasRenderer {
  public:
    struct Rect { int x = 0, y = 0, width = 0, height = 0; };
    // Shared with any UI hit testing. The image itself includes exported axes,
    // a signed colour legend, and explicit numerical sampling budgets.
    static Rect plotRect(int w,int h) {
        const int scale = w>=1050 && h>=650 ? 2 : 1;
        const int left = std::min(w/5,58*scale);
        const int top = std::min(h/3,49*scale);
        const int right = std::min(w/10,18*scale);
        const int bottom = std::min(h/3,73*scale);
        return {left,top,std::max(1,w-left-right),std::max(1,h-top-bottom)};
    }

  private:
    struct Cell { double stretch = 0, area = 0, height = 0; };
    struct RGB { double r = 0, g = 0, b = 0; };
    std::vector<Cell> field_;
    Settings key_;
    uint32_t seed_ = 0;
    bool ready_ = false;
    AtlasStats stats_;

    static bool sameField(const Settings& a,const Settings& b) {
        return a.tilt==b.tilt && a.burnPower==b.burnPower && a.windowPower==b.windowPower &&
               a.kMin==b.kMin && a.kMax==b.kMax && a.tauMin==b.tauMin && a.tauMax==b.tauMax;
    }
    template<class Work> static void parallelRows(int rows,unsigned requestedThreads,Work work) {
        std::atomic<int> next{0};
        std::atomic<bool> failed{false};
        auto run=[&] {
            try { for(int row;!failed.load() && (row=next.fetch_add(1))<rows;) work(row); }
            catch(...) { failed=true; }
        };
        const unsigned available=std::max(1u,std::thread::hardware_concurrency());
        const unsigned wanted=requestedThreads?requestedThreads:std::min(16u,available);
        const unsigned count=std::min(unsigned(rows),std::min(64u,std::max(1u,wanted)));
        std::vector<std::thread> workers;
        workers.reserve(count-1);
        try { for(unsigned i=1;i<count;++i) workers.emplace_back(run); }
        catch(const std::system_error&) { /* Existing workers and caller finish. */ }
        run();
        for(auto& worker:workers) worker.join();
        if(failed) throw std::runtime_error("Spin parameter atlas sampling failed");
    }
    template<class Cancel> bool sample(const Settings& s,uint32_t seed,int cols,int rows,
                                      unsigned threads,Cancel cancelled) {
        if(ready_ && seed_==seed && stats_.gridWidth==cols && stats_.gridHeight==rows && sameField(key_,s))
            return !cancelled();
        const int burn=1<<s.burnPower, window=1<<s.windowPower;
        std::vector<Cell> next(size_t(cols)*size_t(rows));
        const Vec3 start=initialPoint(seed,0), initialTangent=tangentAt(start);
        parallelRows(rows,threads,[&](int y) {
            Settings cellSettings=s;
            cellSettings.damping=s.tauMax-(s.tauMax-s.tauMin)*(double(y)+.5)/rows;
            for(int x=0;x<cols;++x) {
                if(cancelled()) return;
                cellSettings.twist=s.kMin+(s.kMax-s.kMin)*(double(x)+.5)/cols;
                const Model model(cellSettings);
                Vec3 r=start, tangent=initialTangent;
                for(int j=0;j<burn;++j) {
                    if((j&31)==0 && cancelled()) return;
                    const auto step=model.step(r,tangent);
                    r=step.point; tangent=step.tangent;
                }
                double stretch=0, area=0;
                for(int j=0;j<window;++j) {
                    if((j&31)==0 && cancelled()) return;
                    const auto step=model.step(r,tangent);
                    r=step.point; tangent=step.tangent;
                    stretch+=step.logStretch; area+=step.logArea;
                }
                if(!std::isfinite(stretch) || !std::isfinite(area) || !quantum::finite(r))
                    throw std::runtime_error("Nonfinite spin parameter cell");
                next[size_t(y)*cols+x]={stretch/window,area/window,r.z};
            }
        });
        if(cancelled()) return false;
        AtlasStats summary;
        summary.cells=uint64_t(cols)*uint64_t(rows);
        summary.gridWidth=cols; summary.gridHeight=rows; summary.burnIn=burn; summary.window=window;
        summary.minLyapunov=std::numeric_limits<double>::infinity();
        summary.maxLyapunov=-std::numeric_limits<double>::infinity();
        for(const Cell& cell:next) {
            summary.minLyapunov=std::min(summary.minLyapunov,cell.stretch);
            summary.maxLyapunov=std::max(summary.maxLyapunov,cell.stretch);
            summary.meanLyapunov+=cell.stretch; summary.meanAreaRate+=cell.area;
        }
        summary.meanLyapunov/=double(summary.cells); summary.meanAreaRate/=double(summary.cells);
        // Cancelled or failed work never replaces the last complete field.
        if(cancelled()) return false;
        field_.swap(next); key_=s; seed_=seed; stats_=summary; ready_=true;
        return true;
    }
    static RGB hsv(double hue,double saturation,double value) {
        hue-=std::floor(hue);
        const double h=hue*6, f=h-std::floor(h);
        const double p=value*(1-saturation),q=value*(1-saturation*f),t=value*(1-saturation*(1-f));
        switch(int(h)%6) {
        case 0:return {value,t,p}; case 1:return {q,value,p}; case 2:return {p,value,t};
        case 3:return {p,q,value}; case 4:return {t,p,value}; default:return {value,p,q};
        }
    }
    static RGB color(double value,int palette,double range,double hue,double exposure) {
        // A fixed symmetric range per field preserves the visible zero/sign.
        // Exposure changes saturation/value, never sign or the stated range.
        const double u=std::clamp(value/range,-1.0,1.0), strength=std::abs(u);
        if(palette==2) return hsv(hue+.62-.57*(u+1)*.5,.78,.3+.65*std::min(1.0,exposure/1.4));
        // Very small hue modulation keeps positive growth warm and negative
        // growth cool even when the shared palette-phase control changes.
        const double shift=.025*std::sin(hue*2*quantum::pi);
        const double h=(u>=0?.07:.61)+shift;
        return hsv(h,.12+.8*std::pow(strength,.55),
                   .20+.77*std::pow(strength,1.0/std::max(.25,exposure)));
    }
    static uint32_t pack(RGB c) {
        auto channel=[](double v) { return uint32_t(std::clamp(v,0.0,1.0)*255+.5); };
        return 0xff000000u|(channel(c.r)<<16)|(channel(c.g)<<8)|channel(c.b);
    }
    static void pixel(std::vector<uint32_t>& out,int w,int h,int x,int y,uint32_t c) {
        if(x>=0 && y>=0 && x<w && y<h) out[size_t(y)*w+x]=c;
    }
    static void line(std::vector<uint32_t>& out,int w,int h,int x0,int y0,int x1,int y1,uint32_t c) {
        const int steps=std::max(std::abs(x1-x0),std::abs(y1-y0));
        if(steps==0) { pixel(out,w,h,x0,y0,c); return; }
        for(int j=0;j<=steps;++j)
            pixel(out,w,h,int(std::lround(x0+(x1-x0)*double(j)/steps)),
                  int(std::lround(y0+(y1-y0)*double(j)/steps)),c);
    }
    static std::array<uint8_t,7> glyph(char c) {
        switch(c) {
        case 'A':return {14,17,17,31,17,17,17}; case 'B':return {30,17,17,30,17,17,30};
        case 'C':return {14,17,16,16,16,17,14}; case 'D':return {30,17,17,17,17,17,30};
        case 'E':return {31,16,16,30,16,16,31}; case 'F':return {31,16,16,30,16,16,16};
        case 'G':return {14,17,16,23,17,17,15}; case 'H':return {17,17,17,31,17,17,17};
        case 'I':return {14,4,4,4,4,4,14}; case 'J':return {7,2,2,2,2,18,12};
        case 'K':return {17,18,20,24,20,18,17}; case 'L':return {16,16,16,16,16,16,31};
        case 'M':return {17,27,21,21,17,17,17}; case 'N':return {17,25,25,21,19,19,17};
        case 'O':return {14,17,17,17,17,17,14}; case 'P':return {30,17,17,30,16,16,16};
        case 'Q':return {14,17,17,17,21,18,13}; case 'R':return {30,17,17,30,20,18,17};
        case 'S':return {15,16,16,14,1,1,30}; case 'T':return {31,4,4,4,4,4,4};
        case 'U':return {17,17,17,17,17,17,14}; case 'V':return {17,17,17,17,17,10,4};
        case 'W':return {17,17,17,21,21,21,10}; case 'X':return {17,17,10,4,10,17,17};
        case 'Y':return {17,17,10,4,4,4,4}; case 'Z':return {31,1,2,4,8,16,31};
        case '0':return {14,17,19,21,25,17,14}; case '1':return {4,12,4,4,4,4,14};
        case '2':return {14,17,1,2,4,8,31}; case '3':return {30,1,1,14,1,1,30};
        case '4':return {2,6,10,18,31,2,2}; case '5':return {31,16,16,30,1,1,30};
        case '6':return {14,16,16,30,17,17,14}; case '7':return {31,1,2,4,8,8,8};
        case '8':return {14,17,17,14,17,17,14}; case '9':return {14,17,17,15,1,1,14};
        case '-':return {0,0,0,31,0,0,0}; case '+':return {0,4,4,31,4,4,0};
        case '.':return {0,0,0,0,0,12,12}; case ':':return {0,12,12,0,12,12,0};
        case '=':return {0,0,31,0,31,0,0}; case '/':return {1,1,2,4,8,16,16};
        case '(':return {2,4,8,8,8,4,2}; case ')':return {8,4,2,2,2,4,8};
        default:return {};
        }
    }
    static void text(std::vector<uint32_t>& out,int w,int h,int x,int y,const std::string& s,uint32_t c,int scale) {
        for(char ch:s) {
            const auto g=glyph(ch);
            for(int row=0;row<7;++row) for(int col=0;col<5;++col) if(g[row]&(1<<(4-col)))
                for(int dy=0;dy<scale;++dy) for(int dx=0;dx<scale;++dx)
                    pixel(out,w,h,x+col*scale+dx,y+row*scale+dy,c);
            x+=6*scale;
        }
    }
    static std::string number(double value) {
        char buffer[48]; std::snprintf(buffer,sizeof(buffer),"%.3g",value);
        std::string s(buffer); for(char& c:s) if(c=='e') c='E'; return s;
    }

  public:
    AtlasStats stats() const { return stats_; }
    bool render(const Settings& s,const quantum::Settings& q,int palette,double hue,bool guides,
                int w,int h,int quality,unsigned threads,std::vector<uint32_t>& out,
                const std::atomic<uint64_t>* revision=nullptr,uint64_t wantedRevision=0,int historyPreviewPower=0) {
        auto cancelled=[&] { return revision && revision->load()!=wantedRevision; };
        if(cancelled()) return false;
        if(!valid(s) || !quantum::valid(q) || !std::isfinite(hue) || w<1 || h<1)
            throw std::runtime_error("Invalid spin parameter atlas settings");
        palette=std::clamp(palette,0,2);
        const Rect plot=plotRect(w,h);
        int cols=std::min(256,1<<s.gridPower);
        if(quality>1) cols=std::min(512,cols*2);
        if(historyPreviewPower>0) cols=std::min(cols,historyPreviewPower<=17?32:64);
        int rows=std::clamp(int(std::lround(double(cols)*plot.height/plot.width)),16,512);
        // A palette/cursor change usually requests preview, medium and full
        // passes again. Reuse an already denser complete field on every pass;
        // otherwise preview would evict the full cache and trigger needless
        // simulation despite unchanged dynamics. An explicit grid-budget
        // change still starts a field with the newly selected budget.
        if(ready_ && seed_==q.randomSeed && sameField(key_,s) && key_.gridPower==s.gridPower &&
           stats_.gridWidth>=cols && stats_.gridHeight>=rows) {
            cols=stats_.gridWidth; rows=stats_.gridHeight;
        }
        if(!sample(s,q.randomSeed,cols,rows,threads,cancelled)) return false;
        double range=palette==2?1.0:.05;
        if(palette!=2) for(const Cell& cell:field_)
            range=std::max(range,std::abs(palette==0?cell.stretch:cell.area));
        // Nearest-cell reconstruction intentionally avoids inventing smooth
        // dynamics between distinct parameter samples at a basin boundary.
        std::vector<uint32_t> next(size_t(w)*size_t(h),0xff070b14);
        for(int y=0;y<plot.height;++y) {
            if((y&15)==0 && cancelled()) return false;
            const int gy=std::min(rows-1,int(int64_t(y)*rows/plot.height));
            for(int x=0;x<plot.width;++x) {
                const int gx=std::min(cols-1,int(int64_t(x)*cols/plot.width));
                const Cell& cell=field_[size_t(gy)*cols+gx];
                const double value=palette==0?cell.stretch:palette==1?cell.area:cell.height;
                pixel(next,w,h,plot.x+x,plot.y+y,pack(color(value,palette,range,hue,q.exposure)));
            }
        }
        const int scale=w>=1050 && h>=650?2:1, unit=6*scale;
        const uint32_t label=0xffd9e4f2, secondary=0xff8da4bd, border=0xff7991ad;
        if(guides) for(int j=1;j<4;++j) {
            const int xx=plot.x+(plot.width-1)*j/4, yy=plot.y+(plot.height-1)*j/4;
            for(int y=plot.y;y<plot.y+plot.height;y+=4) pixel(next,w,h,xx,y,0xff687385);
            for(int x=plot.x;x<plot.x+plot.width;x+=4) pixel(next,w,h,x,yy,0xff687385);
        }
        line(next,w,h,plot.x-1,plot.y-1,plot.x+plot.width,plot.y-1,border);
        line(next,w,h,plot.x-1,plot.y+plot.height,plot.x+plot.width,plot.y+plot.height,border);
        line(next,w,h,plot.x-1,plot.y-1,plot.x-1,plot.y+plot.height,border);
        line(next,w,h,plot.x+plot.width,plot.y-1,plot.x+plot.width,plot.y+plot.height,border);
        // The selected orbit is a LIVE overlay and is deliberately absent
        // from the field-cache key. Moving it does not resample the atlas.
        if(s.twist>=s.kMin && s.twist<=s.kMax && s.damping>=s.tauMin && s.damping<=s.tauMax) {
            const int mx=plot.x+int(std::lround((s.twist-s.kMin)/(s.kMax-s.kMin)*(plot.width-1)));
            const int my=plot.y+int(std::lround((s.tauMax-s.damping)/(s.tauMax-s.tauMin)*(plot.height-1)));
            const int arm=5*scale;
            for(int thickness=-1;thickness<=1;++thickness) {
                line(next,w,h,mx-arm,my+thickness,mx+arm,my+thickness,0xff000000);
                line(next,w,h,mx+thickness,my-arm,mx+thickness,my+arm,0xff000000);
            }
            line(next,w,h,mx-arm,my,mx+arm,my,0xffffffff);
            line(next,w,h,mx,my-arm,mx,my+arm,0xffffffff);
        }
        const std::string title=palette==0?"FINITE-TIME LYAPUNOV ESTIMATE":palette==1?"MEAN LOG AREA FACTOR":"FINAL HEIGHT Z";
        text(next,w,h,plot.x,7*scale,title,label,scale);
        text(next,w,h,plot.x,20*scale,"ONE INITIAL STATE  TILT="+number(s.tilt)+"  SEED="+std::to_string(q.randomSeed),secondary,scale);
        text(next,w,h,plot.x,33*scale,"GRID="+std::to_string(cols)+"X"+std::to_string(rows)+
             "  BURN="+std::to_string(stats_.burnIn)+"  WINDOW="+std::to_string(stats_.window),secondary,scale);
        text(next,w,h,3*scale,plot.y-10*scale,"TAU",label,scale);
        text(next,w,h,3*scale,plot.y,number(s.tauMax),secondary,scale);
        text(next,w,h,3*scale,plot.y+plot.height/2-3*scale,number((s.tauMin+s.tauMax)*.5),secondary,scale);
        text(next,w,h,3*scale,plot.y+plot.height-7*scale,number(s.tauMin),secondary,scale);
        const int axisY=plot.y+plot.height+7*scale;
        const std::string right=number(s.kMax);
        text(next,w,h,plot.x,axisY,number(s.kMin),secondary,scale);
        text(next,w,h,plot.x+plot.width/2-unit/2,axisY,"K",label,scale);
        text(next,w,h,plot.x+plot.width-int(right.size())*unit,axisY,right,secondary,scale);
        const int barY=axisY+15*scale, barWidth=std::max(1,plot.width), barHeight=6*scale;
        for(int x=0;x<barWidth;++x) {
            const double value=range*(2*(double(x)+.5)/barWidth-1);
            const uint32_t c=pack(color(value,palette,range,hue,q.exposure));
            for(int y=0;y<barHeight;++y) pixel(next,w,h,plot.x+x,barY+y,c);
        }
        const int legendY=barY+barHeight+4*scale;
        text(next,w,h,plot.x,legendY,number(-range),secondary,scale);
        text(next,w,h,plot.x+plot.width/2-unit/2,legendY,"0",secondary,scale);
        const std::string positive=number(range);
        text(next,w,h,plot.x+plot.width-int(positive.size())*unit,legendY,positive,secondary,scale);
        const std::string footer=palette==0?"LOG STRETCH / STEP - FINITE WINDOW":palette==1?"LOG AREA / STEP - FINITE WINDOW":"Z AFTER BURN + WINDOW";
        text(next,w,h,plot.x,legendY+12*scale,footer,secondary,scale);
        if(cancelled()) return false;
        out.swap(next);
        return true;
    }
};
} // namespace spin
