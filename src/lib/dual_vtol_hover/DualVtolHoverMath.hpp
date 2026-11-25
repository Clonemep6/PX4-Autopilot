#pragma once
// Simulator-only dual-rotor hover allocation reference and PX4 adapter kernel.
// Copyright (c) 2026. SPDX-License-Identifier: BSD-3-Clause
// This is a quasi-static, fully-deployed rotor model, NOT a flight-validated model.
#include <array>
#include <cmath>
#include <algorithm>

namespace dvh {
constexpr float pi = 3.14159265358979323846f;
inline float bound(float x, float lo, float hi) { return std::max(lo, std::min(x, hi)); }
using V5 = std::array<float, 5>; // roll/pitch/yaw [Nm], forward/upward [N]
using C5 = std::array<float, 5>; // T_L/Tmax,T_R/Tmax,alpha_L/limit,alpha_R/limit,cyclic
struct Point { float rpm, thrust, torque; };
// Original TOTAL rotor curve; derating applied once in lookup/inverse.
constexpr Point curve[] = {
 {3918,12.50347875f,.25f},{4154,14.17060925f,.29f},{4376,15.79851315f,.32f},
 {4587,17.3969971f,.35f},{4889,19.78001305f,.40f},{5179,22.2022556f,.45f},
 {5464,24.712758f,.50f},{5746,27.33113355f,.55f},{6026,30.09660885f,.61f},
 {6306,33.01899055f,.67f},{6584,36.11789195f,.73f},{6862,39.4031197f,.80f},
 {7138,42.86486715f,.87f},{7412,46.51294095f,.94f},{7685,50.3277278f,1.02f},
 {7956,54.2896144f,1.10f},{8226,58.39860075f,1.19f},{8496,62.61546025f,1.27f},
 {8766,66.9205796f,1.37f},{9035,71.30415215f,1.46f},{9300,75.71714465f,1.55f},
 {10067,89.632781f,1.85f}};
constexpr unsigned ncurve = sizeof(curve) / sizeof(curve[0]);
inline Point lookup(float omega) {
 if (!std::isfinite(omega) || omega <= 0.f) { return {}; }
 const float rpm = omega * 60.f / (2.f*pi);
 if (rpm >= curve[ncurve-1].rpm) { auto p=curve[ncurve-1]; p.thrust*=.9f; p.torque*=.9f; return p; }
 Point a{};
 for (const auto &b : curve) {
  if (rpm <= b.rpm) {
   const float w=(rpm*rpm-a.rpm*a.rpm)/(b.rpm*b.rpm-a.rpm*a.rpm);
   return {rpm,.9f*(a.thrust+w*(b.thrust-a.thrust)),.9f*(a.torque+w*(b.torque-a.torque))};
  }
  a=b;
 }
 return {};
}
inline float omegaForThrust(float thrust) {
 if (!std::isfinite(thrust) || thrust <= 0.f) { return 0.f; }
 const float target=thrust/.9f;
 Point a{};
 for (const auto &b : curve) {
  if (target <= b.thrust) {
   const float w=bound((target-a.thrust)/(b.thrust-a.thrust),0.f,1.f);
   return std::sqrt(a.rpm*a.rpm+w*(b.rpm*b.rpm-a.rpm*a.rpm))*2.f*pi/60.f;
  }
  a=b;
 }
 return curve[ncurve-1].rpm*2.f*pi/60.f;
}
struct Config {
 float omega_max{1000.f}; // MUST match SIM_GZ_EC_MAX1/MAX2; MIN1/MIN2=0
 float arm_l{1.33714025f}, arm_r{1.33695575f};
 float pivot_x{.02577691f}, pivot_z{.01250282f}; // body FRD relative to COM, m
 float cyclic_pitch{-.06194f}, cyclic_roll{.02178f}; // Nm / (N * normalized swash)
 float tilt_limit{10.f*pi/180.f};
 std::array<float,2> tilt_trim{{-.0534036732f,-.0534036732f}}; // joint rad
 std::array<float,3> moment_scale{{10.f,2.f,10.f}}; // Nm per normalized torque input
 float tilt_rate{.3f}; // rad/s command cap; mechanics have their own dynamics
 float cyclic_rate{2.f}; // normalized swash /s
 bool valid() const {
  return std::isfinite(omega_max) && omega_max >= 100.f && omega_max <= 1054.214f
   && std::isfinite(arm_l) && arm_l>.1f && std::isfinite(arm_r) && arm_r>.1f
   && std::isfinite(pivot_x) && std::isfinite(pivot_z)
   && std::isfinite(cyclic_pitch) && std::fabs(cyclic_pitch)>.001f
   && std::isfinite(cyclic_roll) && std::fabs(cyclic_roll)<.2f
   && std::isfinite(tilt_limit) && tilt_limit>.001f && tilt_limit<.524f
   && std::isfinite(tilt_trim[0]) && std::isfinite(tilt_trim[1])
   && std::fabs(tilt_trim[0])+tilt_limit < pi/2.f
   && std::fabs(tilt_trim[1])+tilt_limit < pi/2.f
   && std::isfinite(tilt_rate) && tilt_rate>0.f
   && std::isfinite(cyclic_rate) && cyclic_rate>0.f
   && std::isfinite(moment_scale[0]) && moment_scale[0]>.01f
   && std::isfinite(moment_scale[1]) && moment_scale[1]>.01f
   && std::isfinite(moment_scale[2]) && moment_scale[2]>.01f;
 }
};
inline C5 clipped(C5 x, float minimum_fraction) {
 x[0]=bound(x[0],minimum_fraction,1.f);x[1]=bound(x[1],minimum_fraction,1.f);
 for (unsigned i=2;i<5;i++) { x[i]=bound(x[i],-1.f,1.f); }
 return x;
}
// Averages over a rotor revolution. Cyclic gain/phase require dynamic measurement.
// A nacelle's radial hub offset is parallel to its thrust and adds no r x F term;
// use the pivot x/z, but the hub's y distance, for the force moment arms.
inline V5 forward(const Config &p, const C5 &x) {
 V5 out{};
 const float tmax=lookup(p.omega_max).thrust;
 for (unsigned side=0;side<2;side++) {
  const float T=x[side]*tmax, a=x[2+side]*p.tilt_limit;
  const float h=T*std::sin(a), v=T*std::cos(a);
  const float y=side==0 ? -p.arm_l : p.arm_r;
  const float spin=side==0 ? 1.f : -1.f;
  const float q=lookup(omegaForThrust(T)).torque;
  const float cyc_r=(side==0 ? p.cyclic_roll : -p.cyclic_roll)*T*x[4];
  out[0]+=-y*v+cyc_r*std::cos(a)-spin*q*std::sin(a);
  out[1]+=p.pivot_x*v+p.pivot_z*h+p.cyclic_pitch*T*x[4];
  out[2]+=-y*h+cyc_r*std::sin(a)+spin*q*std::cos(a);
  out[3]+=h;out[4]+=v;
 }
 return out;
}
inline V5 scales(const Config &p) {
 const float f=2.f*lookup(p.omega_max).thrust;
 return {{p.moment_scale[0],p.moment_scale[1],p.moment_scale[2],f,f}};
}
inline float cost(const V5 &a, const V5 &target, const V5 &scale) {
 constexpr float weight[5]={3.f,3.f,1.f,1.f,4.f};
 float sum=0;
 for(unsigned i=0;i<5;i++){const float e=(a[i]-target[i])*weight[i]/scale[i];sum+=e*e;}
 return sum;
}
inline bool solve5(float a[5][5], float b[5], C5 &x) {
 for(unsigned k=0;k<5;k++) {
  unsigned pivot=k;
  for(unsigned i=k+1;i<5;i++) { if(std::fabs(a[i][k])>std::fabs(a[pivot][k])) {pivot=i;} }
  if(!std::isfinite(a[pivot][k]) || std::fabs(a[pivot][k])<1e-12f) {return false;}
  if(pivot!=k){for(unsigned j=k;j<5;j++){std::swap(a[k][j],a[pivot][j]);}std::swap(b[k],b[pivot]);}
  const float divisor=a[k][k];
  for(unsigned j=k;j<5;j++){a[k][j]/=divisor;} b[k]/=divisor;
  for(unsigned i=0;i<5;i++){if(i==k){continue;}const float f=a[i][k];
   for(unsigned j=k;j<5;j++){a[i][j]-=f*a[k][j];}b[i]-=f*b[k];}
 }
 for(unsigned i=0;i<5;i++){x[i]=b[i];if(!std::isfinite(x[i])){return false;}}
 return true;
}
struct Result { C5 controls{}; V5 achieved{}, residual{}; bool valid{false}; };
inline Result allocate(const Config &p,const V5 &target,float idle_omega,const C5 *initial=nullptr) {
 Result out{};
 if(!p.valid() || !std::isfinite(idle_omega)){return out;}
 for(float v:target){if(!std::isfinite(v)){return out;}}
 const float tmax=lookup(p.omega_max).thrust;
 const float minfrac=bound(lookup(idle_omega).thrust/tmax,0.f,1.f);
 C5 x{};
 if(initial){x=*initial;}else{
  x[0]=x[1]=bound(target[4]/(2.f*tmax),minfrac,1.f);
  x[4]=bound((target[1]-p.pivot_x*target[4]-p.pivot_z*target[3])/
   (p.cyclic_pitch*std::max(target[4],1.f)),-1.f,1.f);
 }
 for(float v:x){if(!std::isfinite(v)){return out;}}
 x=clipped(x,minfrac);
 const auto scale=scales(p);
 constexpr float weight[5]={3.f,3.f,1.f,1.f,4.f};
 // Bounded damped Gauss-Newton. Vertical, roll and pitch have higher priority.
 for(unsigned iter=0;iter<18;iter++) {
  const auto f=forward(p,x);const float old_cost=cost(f,target,scale);
  if(old_cost<1e-12f){break;}
  float J[5][5]{};
  for(unsigned col=0;col<5;col++){
   C5 xp=x,xm=x;xp[col]+=.001f;xm[col]-=.001f;
   xp=clipped(xp,minfrac);xm=clipped(xm,minfrac);
   const float d=xp[col]-xm[col];
   const auto fp=forward(p,xp),fm=forward(p,xm);
   for(unsigned row=0;row<5;row++){J[row][col]=d>1e-7f ?
     (fp[row]-fm[row])/d*weight[row]/scale[row] : 0.f;}
  }
  float A[5][5]{},b[5]{};
  for(unsigned i=0;i<5;i++){
   for(unsigned row=0;row<5;row++){
    b[i]+=J[row][i]*(target[row]-f[row])*weight[row]/scale[row];
    for(unsigned j=0;j<5;j++){A[i][j]+=J[row][i]*J[row][j];}
   }
   A[i][i]+=1e-5f;
  }
  C5 step{};if(!solve5(A,b,step)){break;}
  bool improved=false;
  for(float multiplier : {1.f,.5f,.25f,.125f,.0625f,.03125f}){
   auto candidate=x;
   for(unsigned i=0;i<5;i++){candidate[i]+=multiplier*step[i];}
   candidate=clipped(candidate,minfrac);
   if(cost(forward(p,candidate),target,scale)<old_cost){x=candidate;improved=true;break;}
  }
  if(!improved){break;}
 }
 out.controls=x;out.achieved=forward(p,x);out.valid=true;
 for(unsigned i=0;i<5;i++){out.residual[i]=target[i]-out.achieved[i];}
 return out;
}
inline void limitServoRate(const Config &p,float dt,const C5 &previous,C5 &next){
 const float d=bound(dt,0.f,.1f);
 for(unsigned i=2;i<4;i++){next[i]=bound(next[i],previous[i]-p.tilt_rate*d/p.tilt_limit,previous[i]+p.tilt_rate*d/p.tilt_limit);}
 next[4]=bound(next[4],previous[4]-p.cyclic_rate*d,previous[4]+p.cyclic_rate*d);
}
// Position-control force split, NED desired force -> roll plus body x/z forces.
// Desired pitch stays level; lateral acceleration still needs body roll.
struct Split {float roll{0},fx{0},fz{0};};
inline Split splitForce(float north,float east,float down,float yaw){
 if(!std::isfinite(north)||!std::isfinite(east)||!std::isfinite(down)||!std::isfinite(yaw)){return {};}
 const float cy=std::cos(yaw),sy=std::sin(yaw);
 const float fx=cy*north+sy*east,fy=-sy*north+cy*east;
 if(down>=-1e-6f){return {};}
 return {std::atan2(fy,-down),fx,-std::sqrt(fy*fy+down*down)};
}
} // namespace dvh
