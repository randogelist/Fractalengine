#pragma once
#include "spin_math.hpp"
#include "spin_atlas.hpp"
#include <array>
#include <atomic>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <utility>
#include <vector>

namespace spin {
struct Stats {
    uint64_t samples=0;
    int chains=0,burnIn=0,window=0;
    Vec3 mean{0,0,0};
    double meanNorm=0,meanLyapunov=0,meanAreaRate=0;
    bool parameterAtlas=false;
    uint64_t cells=0;
    int gridWidth=0,gridHeight=0,startsPerCell=0;
    double minLyapunov=0,maxLyapunov=0;
};

// Deterministic trajectories: the seed chooses only the 256 initial states.
// Each stored sample carries rolling finite-window tangent/area growth estimates.
// Display choices and atlas bounds do not change the trajectory prefix.
class CloudRenderer {
    static constexpr int chainCount=256;
    struct Point { float x,y,z,growth,area; };
    struct LogPair { double growth=0,area=0; };
    struct Chain {
        Vec3 r{0,0,1},tangent{1,0,0};
        std::vector<LogPair> ring;
        size_t cursor=0;
        double growthSum=0,areaSum=0;
    };
    struct Bin { float count=0,growth=0,area=0; };
    struct RGB { double r=0,g=0,b=0; };
    std::vector<Point> cloud_;
    std::array<Chain,chainCount> chains_{};
    Settings key_;
    uint32_t seed_=0;
    bool ready_=false;
    Stats stats_;
    AtlasRenderer atlas_;

    static bool sameDynamics(const Settings& a,const Settings& b) {
        return a.twist==b.twist && a.tilt==b.tilt && a.damping==b.damping &&
               a.burnPower==b.burnPower && a.windowPower==b.windowPower;
    }
    template<class Work> static void parallelChains(Work work,unsigned requestedThreads) {
        std::atomic<int> next{0};
        std::atomic<bool> failed{false};
        auto run=[&] {
            try { for(int i;!failed && (i=next.fetch_add(1))<chainCount;) work(i); }
            catch(...) { failed=true; }
        };
        const unsigned available=std::max(1u,std::thread::hardware_concurrency());
        const unsigned wanted=requestedThreads?requestedThreads:std::min(16u,available);
        const unsigned count=std::min(64u,std::max(1u,wanted));
        std::vector<std::thread> workers;
        workers.reserve(count-1);
        try { for(unsigned i=1;i<count;++i) workers.emplace_back(run); }
        catch(const std::system_error&) { /* Existing workers and this thread finish. */ }
        run();
        for(auto& worker:workers) worker.join();
        if(failed) throw std::runtime_error("Spin trajectory sampling failed");
    }
    static void advance(Chain& c,const Model& model) {
        const Step step=model.step(c.r,c.tangent);
        c.r=step.point; c.tangent=step.tangent;
        auto& old=c.ring[c.cursor];
        c.growthSum+=step.logStretch-old.growth;
        c.areaSum+=step.logArea-old.area;
        old={step.logStretch,step.logArea};
        c.cursor=(c.cursor+1)%c.ring.size();
        // Re-summing once per ring cycle limits cancellation drift in long runs.
        if(c.cursor==0) {
            c.growthSum=0; c.areaSum=0;
            for(const auto& term:c.ring) { c.growthSum+=term.growth; c.areaSum+=term.area; }
        }
    }
    template<class Cancel> bool sample(const Settings& s,uint32_t seed,size_t wanted,
                                      unsigned threads,Cancel cancelled) {
        const Model model(s);
        const size_t window=size_t(1)<<s.windowPower;
        const int burn=std::max(1<<s.burnPower,1<<s.windowPower);
        if(!ready_ || seed_!=seed || !sameDynamics(key_,s)) {
            // Keep the last complete cache intact until the new cloud finishes.
            std::array<Chain,chainCount> fresh{};
            std::vector<Point> freshCloud(wanted);
            parallelChains([&](int id) {
                Chain c;
                c.r=initialPoint(seed,id); c.tangent=tangentAt(c.r);
                c.ring.resize(window);
                for(int i=0;i<burn;++i) {
                    if((i&63)==0 && cancelled()) return;
                    advance(c,model);
                }
                for(size_t step=0;step<wanted/chainCount;++step) {
                    if((step&63)==0 && cancelled()) return;
                    advance(c,model);
                    freshCloud[step*chainCount+size_t(id)]={float(c.r.x),float(c.r.y),float(c.r.z),
                        float(c.growthSum/window),float(c.areaSum/window)};
                }
                fresh[id]=std::move(c);
            },threads);
            if(cancelled()) return false;
            chains_=std::move(fresh); cloud_=std::move(freshCloud);
            key_=s; seed_=seed; ready_=true;
        }
        if(cloud_.size()>=wanted) return !cancelled();
        const size_t old=cloud_.size(),first=old/chainCount,last=wanted/chainCount;
        const auto saved=chains_;
        cloud_.resize(wanted);
        try {
            parallelChains([&](int id) {
                Chain c=chains_[id];
                for(size_t step=first;step<last;++step) {
                    if((step&63)==0 && cancelled()) return;
                    advance(c,model);
                    cloud_[step*chainCount+size_t(id)]={float(c.r.x),float(c.r.y),float(c.r.z),
                        float(c.growthSum/window),float(c.areaSum/window)};
                }
                chains_[id]=std::move(c);
            },threads);
        } catch(...) { cloud_.resize(old); chains_=saved; throw; }
        if(cancelled()) { cloud_.resize(old); chains_=saved; return false; }
        return true;
    }

    static RGB hsv(double hue,double saturation,double value) {
        hue-=std::floor(hue);
        const double h=hue*6, f=h-std::floor(h);
        const double p=value*(1-saturation), q=value*(1-saturation*f), t=value*(1-saturation*(1-f));
        switch(int(h)%6) {
            case 0:return {value,t,p}; case 1:return {q,value,p}; case 2:return {p,value,t};
            case 3:return {p,q,value}; case 4:return {t,p,value}; default:return {value,p,q};
        }
    }
    static uint32_t pack(RGB c) {
        auto channel=[](double x) { return uint32_t(std::clamp(x,0.0,1.0)*255+.5); };
        return 0xff000000u|(channel(c.r)<<16)|(channel(c.g)<<8)|channel(c.b);
    }
    static void pixel(std::vector<uint32_t>& out,int w,int h,int x,int y,uint32_t color) {
        if(x>=0 && y>=0 && x<w && y<h) out[size_t(y)*w+x]=color;
    }
    static void line(std::vector<uint32_t>& out,int w,int h,int x0,int y0,int x1,int y1,
                     uint32_t color,int clipLeft=0,int clipRight=-1) {
        if(clipRight<0) clipRight=w;
        const int steps=std::max(std::abs(x1-x0),std::abs(y1-y0));
        if(steps==0) { if(x0>=clipLeft && x0<clipRight) pixel(out,w,h,x0,y0,color); return; }
        for(int i=0;i<=steps;++i) {
            const int x=int(std::lround(x0+(x1-x0)*double(i)/steps));
            if(x>=clipLeft && x<clipRight)
                pixel(out,w,h,x,int(std::lround(y0+(y1-y0)*double(i)/steps)),color);
        }
    }
    static std::array<uint8_t,7> glyph(char c) {
        // Small exportable chart labels. Each row is five bits, left to right.
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
        case '-':return {0,0,0,31,0,0,0}; case '.':return {0,0,0,0,0,12,12};
        default:return {};
        }
    }
    static void text(std::vector<uint32_t>& out,int w,int h,int x,int y,const std::string& label,uint32_t color,int scale) {
        for(char c:label) {
            auto g=glyph(c);
            for(int row=0;row<7;++row) for(int col=0;col<5;++col) if(g[row]&(1<<(4-col)))
                for(int j=0;j<scale;++j) for(int i=0;i<scale;++i) pixel(out,w,h,x+col*scale+i,y+row*scale+j,color);
            x+=6*scale;
        }
    }


    static RGB signedColor(double value) {
        const double t=std::clamp(std::abs(value),0.0,1.0);
        const RGB neutral{.85,.86,.87};
        const RGB end=value>=0?RGB{1,.24,.075}:RGB{.12,.48,1};
        return {neutral.r+(end.r-neutral.r)*t,neutral.g+(end.g-neutral.g)*t,
                neutral.b+(end.b-neutral.b)*t};
    }
    static void legend(std::vector<uint32_t>& out,int w,int h,int palette) {
        if(palette==0 || w<260 || h<140) return;
        const int scale=std::max(1,std::min(w/900,h/500));
        const int bw=std::min(240*scale,w-32),left=16,top=h-43*scale;
        for(int y=top-12*scale;y<h-5;++y)
            for(int x=8;x<std::min(w,bw+28);++x) pixel(out,w,h,x,y,0xff070b14);
        text(out,w,h,left,top-10*scale,palette==1?"WINDOW LOG STRETCH":"WINDOW LOG AREA",0xffdce4ef,scale);
        for(int x=0;x<bw;++x) for(int y=0;y<8*scale;++y)
            pixel(out,w,h,left+x,top+y,pack(signedColor(2.0*x/(bw-1)-1)));
        text(out,w,h,left,top+11*scale,"-1",0xffb8c4d6,scale);
        text(out,w,h,left+bw/2-3*scale,top+11*scale,"0",0xffb8c4d6,scale);
        text(out,w,h,left+bw-6*scale,top+11*scale,"1",0xffb8c4d6,scale);
        text(out,w,h,left,top+22*scale,"NATS PER STEP - CLIPPED",0xffb8c4d6,scale);
    }

  public:
    Stats stats() const { return stats_; }
    bool render(int generation,const Settings& s,const quantum::Settings& q,int plane,int palette,
                double hue,bool guides,int w,int h,int quality,unsigned threads,
                std::vector<uint32_t>& out,const std::atomic<uint64_t>* revision=nullptr,
                uint64_t wantedRevision=0,int historyPreviewPower=0) {
        auto cancelled=[&] { return revision && revision->load()!=wantedRevision; };
        if(cancelled()) return false;
        if(!valid(s) || !quantum::valid(q) || !std::isfinite(hue) || generation<1 || generation>3 ||
           plane<0 || plane>2 || palette<0 || palette>2 || w<1 || h<1)
            throw std::runtime_error("Invalid spin render settings");
        if(generation==3 && plane==0) {
            if(!atlas_.render(s,q,palette,hue,guides,w,h,quality,threads,out,revision,wantedRevision,historyPreviewPower)) return false;
            const auto a=atlas_.stats();
            Stats next;
            next.parameterAtlas=true; next.cells=a.cells; next.gridWidth=a.gridWidth; next.gridHeight=a.gridHeight;
            next.burnIn=a.burnIn; next.window=a.window; next.startsPerCell=a.startsPerCell;
            next.minLyapunov=a.minLyapunov; next.maxLyapunov=a.maxLyapunov;
            next.meanLyapunov=a.meanLyapunov; next.meanAreaRate=a.meanAreaRate;
            stats_=next;
            return true;
        }
        size_t count=size_t(1)<<q.samplePower;
        if(quality>1) count=std::min(size_t(1)<<23,count*4);
        if(historyPreviewPower>0) count=std::min(count,size_t(1)<<std::clamp(historyPreviewPower,14,23));
        if(!sample(s,q.randomSeed,count,threads,cancelled)) return false;
        // The third entry reserves plane 0 for the parameter atlas.
        const int view=generation==3?plane-1:plane;
        const bool atlas=view==1,pair=view==2;
        const int split=pair?w/2:w,panels=pair?2:1;
        const double cy=(h-1)*.5;
        const double cosineY=std::cos(q.yaw),sineY=std::sin(q.yaw),cosineP=std::cos(q.pitch),sineP=std::sin(q.pitch);
        auto camera=[&](Vec3 v) {
            const double x=cosineY*v.x+sineY*v.z,z=-sineY*v.x+cosineY*v.z;
            return Vec3{x,cosineP*v.y-sineP*z,sineP*v.y+cosineP*z};
        };
        auto panelLeft=[&](int panel) { return panel==0?0:split; };
        auto panelWidth=[&](int panel) { return panel==0?split:w-split; };
        auto radius=[&](int panel) { return .46*std::min(panelWidth(panel),h)*q.zoom; };
        auto center=[&](int panel) { return panelLeft(panel)+(panelWidth(panel)-1)*.5; };
        auto projectCamera=[&](Vec3 v,int panel,double& x,double& y) {
            if(atlas) {
                x=(w-1)*.5+std::atan2(v.x,v.z)/(2*quantum::pi)*w*q.zoom;
                y=cy-v.y*.5*h*q.zoom;
                return true;
            }
            x=center(panel)+(panel==0?v.x:-v.x)*radius(panel);
            y=cy-v.y*radius(panel);
            return panel==0?v.z>=0:v.z<=0;
        };
        std::vector<Bin> bins(size_t(w)*h);
        Stats next;
        next.samples=count; next.chains=chainCount;
        next.burnIn=std::max(1<<s.burnPower,1<<s.windowPower); next.window=1<<s.windowPower;
        next.minLyapunov=std::numeric_limits<double>::infinity();
        next.maxLyapunov=-std::numeric_limits<double>::infinity();
        Vec3 sum{0,0,0}; double growthSum=0,areaSum=0;
        for(size_t i=0;i<count;++i) {
            if((i&4095)==0 && cancelled()) return false;
            const Point& p=cloud_[i];
            const Vec3 v{p.x,p.y,p.z};
            sum=sum+v; growthSum+=p.growth; areaSum+=p.area;
            next.minLyapunov=std::min(next.minLyapunov,double(p.growth));
            next.maxLyapunov=std::max(next.maxLyapunov,double(p.growth));
            const Vec3 rotated=camera(v);
            const int panel=pair && rotated.z<0?1:0;
            double x,y;
            if(!projectCamera(rotated,panel,x,y) || x<-1 || y<-1 || x>w || y>h) continue;
            const int ix=int(std::floor(x)),iy=int(std::floor(y));
            const double fx=x-ix,fy=y-iy;
            for(int sy=0;sy<2;++sy) for(int sx=0;sx<2;++sx) {
                int px=ix+sx;
                const int py=iy+sy;
                if(atlas && q.zoom==1) px=(px%w+w)%w;
                if(px<panelLeft(panel) || px>=panelLeft(panel)+panelWidth(panel) || py<0 || py>=h) continue;
                const float weight=float((sx?fx:1-fx)*(sy?fy:1-fy));
                auto& bin=bins[size_t(py)*w+px];
                bin.count+=weight; bin.growth+=weight*p.growth; bin.area+=weight*p.area;
            }
        }
        out.assign(size_t(w)*h,0xff070b14);
        const double logScale=std::log1p(q.exposure*30);
        for(int y=0;y<h;++y) {
            if((y&31)==0 && cancelled()) return false;
            for(int x=0;x<w;++x) {
                const int panel=pair && x>=split?1:0;
                const double rad=radius(panel);
                if(!(rad>0)) continue;
                const auto& bin=bins[size_t(y)*w+x];
                const double xx=(x-center(panel))/rad,yy=(y-cy)/rad,rr=xx*xx+yy*yy;
                if(!atlas && rr>1.012) continue;
                RGB color=atlas?RGB{.019,.029,.052}:RGB{.025+.018*std::sqrt(std::max(0.0,1-rr)),.037,.065};
                const double reference=std::max(.01,atlas?double(count)/(w*double(h)*q.zoom*q.zoom):double(count)/(2*quantum::pi*rad*rad));
                if(bin.count>0) {
                    const double t=std::clamp(std::log1p(q.exposure*bin.count/reference)/logScale,0.0,1.0);
                    const double brightness=.08+.92*std::pow(t,.72);
                    if(palette==0) color=hsv(hue+.73-.69*t,.88-.53*t,brightness);
                    else {
                        color=signedColor((palette==1?bin.growth:bin.area)/bin.count);
                        color={color.r*brightness,color.g*brightness,color.b*brightness};
                    }
                }
                out[size_t(y)*w+x]=pack(color);
            }
        }
        if(guides) {
            auto drawCurve=[&](auto curve) {
                for(int panel=0;panel<panels;++panel) {
                    bool prior=false; double oldX=0,oldY=0;
                    for(int step=0;step<=480;++step) {
                        double x,y;
                        const bool visible=projectCamera(camera(curve(step*(2*quantum::pi/480))),panel,x,y);
                        if(visible && prior && std::abs(x-oldX)<panelWidth(panel)*.25 && std::abs(y-oldY)<h*.25)
                            line(out,w,h,int(std::lround(oldX)),int(std::lround(oldY)),int(std::lround(x)),int(std::lround(y)),
                                 0xff334459,panelLeft(panel),panelLeft(panel)+panelWidth(panel));
                        prior=visible; oldX=x; oldY=y;
                    }
                }
            };
            for(int k=-2;k<=2;++k) {
                const double latitude=k*quantum::pi/6,c=std::cos(latitude),sn=std::sin(latitude);
                drawCurve([=](double t) { return Vec3{c*std::sin(t),sn,c*std::cos(t)}; });
            }
            for(int k=0;k<6;++k) {
                const double longitude=k*quantum::pi/6;
                drawCurve([=](double t) { return Vec3{std::sin(t)*std::sin(longitude),std::cos(t),std::sin(t)*std::cos(longitude)}; });
            }
        }
        if(pair && w>=240 && h>=80) {
            const int scale=std::max(1,std::min(w/720,h/360));
            text(out,w,h,12,12,"FRONT",0xffcbd5e1,scale);
            text(out,w,h,split+12,12,"BACK",0xffcbd5e1,scale);
            line(out,w,h,split,0,split,h-1,0xff243244);
        }
        legend(out,w,h,palette);
        if(cancelled()) return false;
        next.mean=sum/double(count); next.meanNorm=quantum::length(next.mean);
        next.meanLyapunov=growthSum/count; next.meanAreaRate=areaSum/count;
        stats_=next;
        return true;
    }
};
} // namespace spin
