#pragma once
#include "quantum_math.hpp"
#include <array>
#include <atomic>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <vector>

namespace quantum {
struct Stats {
    uint64_t samples = 0;
    int chains = 0, burnIn = 0;
    Vec3 mean{0,0,0};
    double meanNorm = 0;
};

// Sampling is independent of pixels, camera, palette and worker count. Each
// chain has its own reproducible stream. Larger budgets extend the same cloud
// prefix; only the dynamics controls or random seed discard the cache.
class CloudRenderer {
    static constexpr int chainCount = 256, burnCount = 256, meanSteps = 128;
    struct Point { float x, y, z; uint8_t outcome; uint8_t pad[3]{}; };
    struct Chain { Vec3 r{0,0,1}; uint64_t random = 0; };
    struct Bin { float count=0, r=0, g=0, b=0; };
    struct RGB { double r=0, g=0, b=0; };
    std::vector<Point> cloud_;
    std::array<Chain,chainCount> chains_{};
    std::array<Vec3,meanSteps+1> transient_{};
    Settings key_;
    bool ready_ = false;
    Stats stats_;

    static uint64_t random64(uint64_t& state) {
        uint64_t z=(state+=UINT64_C(0x9e3779b97f4a7c15));
        z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
        z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
        return z^(z>>31);
    }
    static bool sameDynamics(const Settings& a,const Settings& b) {
        return a.strength==b.strength && a.kick==b.kick &&
            a.axisAzimuth==b.axisAzimuth && a.axisElevation==b.axisElevation &&
            a.geometry==b.geometry && a.randomSeed==b.randomSeed;
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
        std::vector<std::thread> threads;
        threads.reserve(count-1);
        try { for(unsigned i=1;i<count;++i) threads.emplace_back(run); }
        catch(const std::system_error&) { /* Existing workers and this thread finish. */ }
        run();
        for(auto& t:threads) t.join();
        if(failed) throw std::runtime_error("Quantum trajectory sampling failed");
    }
    // Zero-probability outcomes are explicitly skipped, including the
    // antipodal projective outcome at strength=1.
    static uint8_t advance(Chain& chain,const Settings& q,const std::vector<Vec3>& dirs,Vec3 axis) {
        Vec3 r=rotate(chain.r,axis,q.kick);
        if(!finite(r)) throw std::runtime_error("Nonfinite quantum state");
        std::array<double,20> probabilities{};
        double total=0;
        int last=-1;
        for(size_t i=0;i<dirs.size();++i) {
            const double p=probability(r,dirs[i],q.strength,int(dirs.size()));
            if(!std::isfinite(p) || p<0) throw std::runtime_error("Invalid Born probability");
            probabilities[i]=p;
            total+=p;
            if(p>0) last=int(i);
        }
        if(!(total>0) || !std::isfinite(total) || last<0)
            throw std::runtime_error("Invalid Born normalization");
        const double target=double(random64(chain.random)>>11)*(1.0/9007199254740992.0)*total;
        double cumulative=0;
        int selected=last;
        for(size_t i=0;i<dirs.size();++i) {
            if(probabilities[i]<=0) continue;
            cumulative+=probabilities[i];
            if(target<cumulative) { selected=int(i); break; }
        }
        chain.r=measure(r,dirs[size_t(selected)],q.strength);
        const double norm2=dot(chain.r,chain.r);
        if(!finite(chain.r) || !std::isfinite(norm2) || !(norm2>0))
            throw std::runtime_error("Invalid conditional quantum state");
        // In the zero-strength rotation limit measure is exactly the identity;
        // correct accumulated floating-point drift without adding noise.
        if(std::abs(norm2-1)>1e-12) chain.r=chain.r/std::sqrt(norm2);
        return uint8_t(selected);
    }
    template<class Cancel> bool sample(const Settings& q,size_t wanted,unsigned threads,Cancel cancelled) {
        const auto& dirs=directions(q.geometry);
        const Vec3 axis=kickAxis(q);
        if(!ready_ || !sameDynamics(key_,q)) {
            ready_=false;
            cloud_.clear();
            std::vector<Vec3> trace(size_t(chainCount)*(meanSteps+1),Vec3{0,0,0});
            parallelChains([&](int id) {
                Chain c;
                c.random=(uint64_t(q.randomSeed)<<32)^UINT64_C(0xd1b54a32d192ed03)^uint64_t(id);
                random64(c.random);
                trace[size_t(id)*(meanSteps+1)]={0,0,1};
                for(int step=1;step<=burnCount;++step) {
                    if((step&31)==0 && cancelled()) return;
                    advance(c,q,dirs,axis);
                    if(step<=meanSteps) trace[size_t(id)*(meanSteps+1)+step]=c.r;
                }
                chains_[id]=c;
            },threads);
            if(cancelled()) return false;
            transient_.fill(Vec3{0,0,0});
            for(int id=0;id<chainCount;++id)
                for(int step=0;step<=meanSteps;++step)
                    transient_[step]=transient_[step]+trace[size_t(id)*(meanSteps+1)+step]/double(chainCount);
            key_=q;
            ready_=true;
        }
        if(cloud_.size()>=wanted) return !cancelled();
        const size_t old=cloud_.size(), firstStep=old/chainCount, lastStep=wanted/chainCount;
        const auto saved=chains_;
        cloud_.resize(wanted);
        try {
            parallelChains([&](int id) {
                Chain c=chains_[id];
                for(size_t step=firstStep;step<lastStep;++step) {
                    if((step&63)==0 && cancelled()) return;
                    const uint8_t outcome=advance(c,q,dirs,axis);
                    cloud_[step*chainCount+size_t(id)]={float(c.r.x),float(c.r.y),float(c.r.z),outcome,{0,0,0}};
                }
                chains_[id]=c;
            },threads);
        } catch(...) {
            cloud_.resize(old); chains_=saved; throw;
        }
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
    static void line(std::vector<uint32_t>& out,int w,int h,int x0,int y0,int x1,int y1,uint32_t color,int clipRight=-1) {
        // Bound work even if an off-screen projected guide is very large.
        const int dx=std::abs(x1-x0),dy=std::abs(y1-y0),steps=std::max(dx,dy);
        if(clipRight<0) clipRight=w;
        if(steps==0) { if(x0<clipRight) pixel(out,w,h,x0,y0,color); return; }
        for(int i=0;i<=steps;++i) {
            const int x=int(std::lround(x0+(x1-x0)*double(i)/steps));
            if(x<clipRight) pixel(out,w,h,x,int(std::lround(y0+(y1-y0)*double(i)/steps)),color);
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
    void meanChart(const Settings& q,int w,int h,int split,std::vector<uint32_t>& out) const {
        const int panel=w-split, scale=std::max(1,std::min(w/620,h/260));
        const int left=split+std::max(25,panel/12), right=w-std::max(15,panel/20);
        const int top=std::max(38,h/5), bottom=h-std::max(58,h/5);
        if(right<=left || bottom<=top) return;
        const uint32_t white=0xffcbd5e1, grid=0xff243244, exact=0xff58e1ed, empirical=0xffffb35b;
        text(out,w,h,left,std::max(10,top-28*scale),"ENSEMBLE MEAN",white,scale);
        for(int i=0;i<=4;++i) {
            const int y=bottom-(bottom-top)*i/4;
            line(out,w,h,left,y,right,y,grid);
            const int x=left+(right-left)*i/4;
            line(out,w,h,x,top,x,bottom,grid);
        }
        text(out,w,h,left-12*scale,top-3*scale,"1",white,scale);
        text(out,w,h,left-12*scale,bottom-3*scale,"0",white,scale);
        text(out,w,h,left,bottom+8*scale,"0",white,scale);
        text(out,w,h,right-18*scale,bottom+8*scale,"128",white,scale);
        double expected=1;
        int px=left, py=top, ex=left, ey=top;
        for(int n=1;n<=meanSteps;++n) {
            const int x=left+int((right-left)*double(n)/meanSteps);
            expected*=meanFactor(q.strength);
            const int y=bottom-int(std::clamp(length(transient_[n]),0.0,1.0)*(bottom-top));
            const int e=bottom-int(std::clamp(expected,0.0,1.0)*(bottom-top));
            line(out,w,h,px,py,x,y,empirical); line(out,w,h,ex,ey,x,e,exact);
            px=ex=x; py=y; ey=e;
        }
        const int legendY=bottom+24*scale;
        text(out,w,h,left,legendY,"EXACT MEAN NORM",exact,scale);
        text(out,w,h,left,legendY+12*scale,"256 CHAINS",empirical,scale);
        if(panel>210*scale)
            text(out,w,h,left,legendY+24*scale,"UPDATES FROM NORTH",white,scale);
    }

  public:
    Stats stats() const { return stats_; }
    bool render(const Settings& q,int plane,int palette,double hue,bool guides,int w,int h,int quality,
                unsigned threads,std::vector<uint32_t>& out,const std::atomic<uint64_t>* revision=nullptr,uint64_t wantedRevision=0,
                int historyPreviewPower=0) {
        auto cancelled=[&] { return revision && revision->load()!=wantedRevision; };
        if(cancelled()) return false;
        if(!valid(q) || w<1 || h<1) throw std::runtime_error("Invalid quantum render settings");
        size_t count=size_t(1)<<q.samplePower;
        if(quality>1) count=std::min(size_t(1)<<23,count*4);
        // Only explicitly provisional frames use a reduced history budget.
        // Small final windows and command-line exports honor the chosen budget.
        if(historyPreviewPower>0) count=std::min(count,size_t(1)<<std::clamp(historyPreviewPower,14,23));
        if(!sample(q,count,threads,cancelled)) return false;
        const bool atlas=plane==1, comparison=plane==2;
        const int projectionW=comparison?std::max(1,int(w*.56)):w;
        const double cx=(projectionW-1)*.5,cy=(h-1)*.5;
        const double radius=.46*std::min(projectionW,h)*q.zoom;
        const double cosineY=std::cos(q.yaw),sineY=std::sin(q.yaw),cosineP=std::cos(q.pitch),sineP=std::sin(q.pitch);
        auto camera=[&](Vec3 v) {
            const double x=cosineY*v.x+sineY*v.z, z=-sineY*v.x+cosineY*v.z;
            return Vec3{x,cosineP*v.y-sineP*z,sineP*v.y+cosineP*z};
        };
        auto project=[&](Vec3 v,double& x,double& y) {
            v=camera(v);
            if(atlas) {
                x=cx+std::atan2(v.x,v.z)/(2*pi)*projectionW*q.zoom;
                y=cy-v.y*.5*h*q.zoom;
                return true;
            }
            x=cx+v.x*radius; y=cy-v.y*radius;
            return v.z>=0;
        };
        std::vector<Bin> bins(size_t(projectionW)*h);
        std::array<RGB,20> outcomeColors{};
        const auto& dirs=directions(q.geometry);
        for(size_t i=0;i<dirs.size();++i) outcomeColors[i]=hsv(hue+double(i)/dirs.size(),.7,1);
        Vec3 sum{0,0,0};
        for(size_t i=0;i<count;++i) {
            if((i&4095)==0 && cancelled()) return false;
            const Point& p=cloud_[i];
            const Vec3 v{p.x,p.y,p.z};
            sum=sum+v;
            double x,y;
            if(!project(v,x,y) || x<-1 || y<-1 || x>projectionW || y>h) continue;
            const int ix=int(std::floor(x)),iy=int(std::floor(y));
            const double fx=x-ix,fy=y-iy;
            const RGB color=palette==1?outcomeColors[p.outcome]:palette==2?hsv(hue+.65-(v.z+1)*.29,.75,1):RGB{};
            for(int sy=0;sy<2;++sy) for(int sx=0;sx<2;++sx) {
                int px=ix+sx;
                const int py=iy+sy;
                // A full atlas is periodic in longitude. Wrap the two edge
                // splats instead of losing half their weight at the seam.
                if(atlas && q.zoom==1) px=(px%projectionW+projectionW)%projectionW;
                if(px<0 || py<0 || px>=projectionW || py>=h) continue;
                const float weight=float((sx?fx:1-fx)*(sy?fy:1-fy));
                auto& bin=bins[size_t(py)*projectionW+px];
                bin.count+=weight;
                if(palette) { bin.r+=weight*float(color.r);bin.g+=weight*float(color.g);bin.b+=weight*float(color.b); }
            }
        }
        out.assign(size_t(w)*h,0xff070b14);
        const double reference=std::max(.01,atlas?double(count)/(projectionW*double(h)*q.zoom*q.zoom):double(count)/(2*pi*radius*radius));
        const double logScale=std::log1p(q.exposure*30);
        for(int y=0;y<h;++y) {
            if((y&31)==0 && cancelled()) return false;
            for(int x=0;x<projectionW;++x) {
                const auto& bin=bins[size_t(y)*projectionW+x];
                const double xx=(x-cx)/radius,yy=(y-cy)/radius,rr=xx*xx+yy*yy;
                if(!atlas && rr>1.012) continue;
                RGB color=atlas?RGB{.019,.029,.052}:RGB{.025+.018*std::sqrt(std::max(0.0,1-rr)),.037,.065};
                if(bin.count>0) {
                    const double t=std::clamp(std::log1p(q.exposure*bin.count/reference)/logScale,0.0,1.0);
                    const double brightness=.08+.92*std::pow(t,.72);
                    if(palette==0) color=hsv(hue+.73-.69*t,.88-.53*t,brightness);
                    else color={bin.r/bin.count*brightness,bin.g/bin.count*brightness,bin.b/bin.count*brightness};
                }
                out[size_t(y)*w+x]=pack(color);
            }
        }
        if(guides) {
            const uint32_t grid=0xff334459;
            auto drawCurve=[&](auto curve) {
                bool prior=false;double oldX=0,oldY=0;
                for(int step=0;step<=240;++step) {
                    double x,y; const bool visible=project(curve(step*(2*pi/240)),x,y);
                    if(visible && prior && std::abs(x-oldX)<projectionW*.25 && std::abs(y-oldY)<h*.25)
                        line(out,w,h,int(std::lround(oldX)),int(std::lround(oldY)),int(std::lround(x)),int(std::lround(y)),grid,projectionW);
                    prior=visible;oldX=x;oldY=y;
                }
            };
            for(int k=-2;k<=2;++k) {
                const double latitude=k*pi/6,c=std::cos(latitude),s=std::sin(latitude);
                drawCurve([=](double t) { return Vec3{c*std::sin(t),s,c*std::cos(t)}; });
            }
            for(int k=0;k<6;++k) {
                const double longitude=k*pi/6;
                drawCurve([=](double t) { return Vec3{std::sin(t)*std::sin(longitude),std::cos(t),std::sin(t)*std::cos(longitude)}; });
            }
            for(Vec3 direction:dirs) {
                double x,y;if(!project(direction,x,y)) continue;
                const int px=int(std::lround(x)),py=int(std::lround(y)),size=std::max(3,w/300);
                line(out,w,h,px-size,py,px+size,py,0xffe3edf7,projectionW);
                line(out,w,h,px,py-size,px,py+size,0xffe3edf7,projectionW);
            }
        }
        if(comparison) meanChart(q,w,h,projectionW,out);
        if(cancelled()) return false;
        Stats next;
        next.samples=count;next.chains=chainCount;next.burnIn=burnCount;
        next.mean=sum/double(count);next.meanNorm=length(next.mean);
        stats_=next;
        return true;
    }
};
} // namespace quantum
