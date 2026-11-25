#include "ActuatorEffectivenessDualVtol.hpp"
#include <px4_platform_common/log.h>
#include <cstdio>
#include <cmath>
using matrix::Vector3f;

ActuatorEffectivenessDualVtol::ActuatorEffectivenessDualVtol(ModuleParams *parent) : ModuleParams(parent)
{
 updateParams();
 PX4_WARN("Dual VTOL hover: experimental SITL allocator; no deployment feedback");
}
float ActuatorEffectivenessDualVtol::readFloat(const char *name, float fallback)
{
 const param_t h=param_find(name);
 if(h==PARAM_INVALID || param_get(h,&fallback)!=0){_config_valid=false;}
 return fallback;
}
int32_t ActuatorEffectivenessDualVtol::readInt(const char *name, int32_t fallback)
{
 const param_t h=param_find(name);
 if(h==PARAM_INVALID || param_get(h,&fallback)!=0){_config_valid=false;}
 return fallback;
}
void ActuatorEffectivenessDualVtol::updateParams()
{
 ModuleParams::updateParams();
 _config_valid=true;
 _run_requested=readInt("DVH_RUN",0)==1;
 _cfg.omega_max=readFloat("DVH_OMEGA_MAX",1000.f);
 _cfg.arm_l=readFloat("DVH_ARM_L",1.33714025f);
 _cfg.arm_r=readFloat("DVH_ARM_R",1.33695575f);
 _cfg.pivot_x=readFloat("DVH_PIV_X",.02577691f);
 _cfg.pivot_z=readFloat("DVH_PIV_Z",.01250282f);
 _cfg.cyclic_pitch=readFloat("DVH_CY_P",-.06194f);
 _cfg.cyclic_roll=readFloat("DVH_CY_R",.02178f);
 _cfg.tilt_limit=readFloat("DVH_TILT_MAX",10.f)*dvh::pi/180.f;
 _cfg.tilt_trim[0]=readFloat("DVH_TRIM_L",-.0534036732f);
 _cfg.tilt_trim[1]=readFloat("DVH_TRIM_R",-.0534036732f);
 _cfg.tilt_rate=readFloat("DVH_TILT_RATE",.3f);
 _cfg.cyclic_rate=readFloat("DVH_CY_RATE",2.f);
 _cfg.moment_scale[0]=readFloat("DVH_MX_SCALE",10.f);
 _cfg.moment_scale[1]=readFloat("DVH_MY_SCALE",2.f);
 _cfg.moment_scale[2]=readFloat("DVH_MZ_SCALE",10.f);
 _idle_omega=readFloat("DVH_IDLE_W",400.f);
 _spool_seconds=readFloat("COM_SPOOLUP_TIME",5.f);
 _config_valid=_config_valid && _cfg.valid() && std::isfinite(_idle_omega)
  && _idle_omega>=0.f && _idle_omega<_cfg.omega_max
  && std::isfinite(_spool_seconds) && _spool_seconds>=1.f;
 // The allocator already performs the thrust -> RPM inverse. Do not transform twice.
 _config_valid=(std::fabs(readFloat("THR_MDL_FAC",0.f))<1e-6f) && _config_valid;
 _config_valid=(readInt("CA_R_REV",0)==0) && _config_valid;
 const float expected_angle[8]={.45f*180.f/dvh::pi,.45f*180.f/dvh::pi,80.f,
   .0314f*180.f/dvh::pi,80.f,90.f,.0314f*180.f/dvh::pi,90.f};
 char name[32];
 for(int i=1;i<=8;i++){
  std::snprintf(name,sizeof(name),"SIM_GZ_SV_FUNC%d",i);
  _config_valid=(readInt(name,0)==200+i) && _config_valid;
  std::snprintf(name,sizeof(name),"SIM_GZ_SV_MINA%d",i);
  _config_valid=(std::fabs(readFloat(name,0.f)+expected_angle[i-1])<.002f) && _config_valid;
  std::snprintf(name,sizeof(name),"SIM_GZ_SV_MAXA%d",i);
  _config_valid=(std::fabs(readFloat(name,0.f)-expected_angle[i-1])<.002f) && _config_valid;
  std::snprintf(name,sizeof(name),"CA_SV%d_SLEW",i-1);
  _config_valid=(std::fabs(readFloat(name,0.f))<1e-6f) && _config_valid;
 }
 for(int i=1;i<=3;i++){
  std::snprintf(name,sizeof(name),"SIM_GZ_EC_FUNC%d",i);
  _config_valid=(readInt(name,0)==100+i) && _config_valid;
  std::snprintf(name,sizeof(name),"SIM_GZ_EC_MIN%d",i);
  _config_valid=(readInt(name,0)==0) && _config_valid;
  std::snprintf(name,sizeof(name),"CA_R%d_SLEW",i-1);
  _config_valid=(std::fabs(readFloat(name,0.f))<1e-6f) && _config_valid;
  if(i<=2){
   std::snprintf(name,sizeof(name),"SIM_GZ_EC_MAX%d",i);
   _config_valid=(std::fabs(float(readInt(name,0))-_cfg.omega_max)<.1f) && _config_valid;
  }
 }
 // This prototype uses the MC stack; no VTOL transition or fixed-wing mode is implemented.
 if(!_config_valid && !_reported_invalid){
  PX4_ERR("DVH config mismatch: outputs inhibited; check airframe wiring/scales and DVH params");
  _reported_invalid=true;
 }else if(_config_valid && _reported_invalid){PX4_INFO("DVH config valid");_reported_invalid=false;}
}

bool ActuatorEffectivenessDualVtol::getEffectivenessMatrix(Configuration &c, EffectivenessUpdateReason reason)
{
 if(reason==EffectivenessUpdateReason::NO_EXTERNAL_UPDATE){return false;}
 // A local Jacobian is supplied for topology/status display. Actual bounded nonlinear
 // allocation is performed in updateSetpoint, as with PX4's helicopter allocators.
 const dvh::Config cfg=_cfg.valid() ? _cfg : dvh::Config{};
 const auto scale=dvh::scales(cfg);
 const float tmax=dvh::lookup(cfg.omega_max).thrust;
 dvh::C5 x{{.364f,.364f,0.f,0.f,.416f}};
 dvh::V5 deriv[5]{};
 for(int col=0;col<5;col++){
  auto a=x,b=x;a[col]+=.001f;b[col]-=.001f;
  const auto fa=dvh::forward(cfg,a),fb=dvh::forward(cfg,b);
  for(int row=0;row<5;row++){deriv[col][row]=(fa[row]-fb[row])/.002f/scale[row];}
 }
 const float w=dvh::omegaForThrust(x[0]*tmax);
 const float dtdw=(dvh::lookup(w+1.f).thrust-dvh::lookup(w-1.f).thrust)*.5f;
 auto add=[&](ActuatorType type,int col,float factor){
  dvh::V5 v{};if(col>=0){for(int row=0;row<5;row++){v[row]=deriv[col][row]*factor;}}
  return c.addActuator(type,Vector3f{v[0],v[1],v[2]},Vector3f{v[3],0.f,-v[4]})>=0;
 };
 bool ok=add(ActuatorType::MOTORS,0,dtdw*cfg.omega_max/tmax);
 ok=add(ActuatorType::MOTORS,1,dtdw*cfg.omega_max/tmax) && ok;
 ok=add(ActuatorType::MOTORS,-1,0.f) && ok; // rear pusher: present but off in hover
 for(int i=0;i<8;i++){
  int col=-1;float factor=0.f;
  if(i==3 || i==6){col=4;factor=.5f;}
  if(i==5 || i==7){col=i==5 ? 2 : 3;factor=(dvh::pi/2.f)/cfg.tilt_limit;}
  ok=add(ActuatorType::SERVOS,col,factor) && ok;
 }
 return ok;
}
void ActuatorEffectivenessDualVtol::setSafeOutput(ActuatorVector &out) const
{
 out.setZero();
 out(3+5)=std::isfinite(_cfg.tilt_trim[0]) ? dvh::bound(_cfg.tilt_trim[0]/(dvh::pi/2.f),-1.f,1.f) : 0.f;
 out(3+7)=std::isfinite(_cfg.tilt_trim[1]) ? dvh::bound(_cfg.tilt_trim[1]/(dvh::pi/2.f),-1.f,1.f) : 0.f;
}
void ActuatorEffectivenessDualVtol::putResidual(const matrix::Vector<float,NUM_AXES> &request,const dvh::V5 &achieved)
{
 const auto scale=dvh::scales(_cfg.valid() ? _cfg : dvh::Config{});
 for(int i=0;i<3;i++){_residual[i]=request(i)-achieved[i]/scale[i];}
 _residual[3]=request(3)-achieved[3]/scale[3];
 _residual[4]=request(4); // no independent body lateral force actuator
 _residual[5]=request(5)+achieved[4]/scale[4];
 for(float &v:_residual){if(!std::isfinite(v)){v=0.f;}else if(std::fabs(v)<1e-5f){v=0.f;}}
}
void ActuatorEffectivenessDualVtol::updateSetpoint(const matrix::Vector<float,NUM_AXES> &request,int matrix_index,
  ActuatorVector &out,const ActuatorVector &minimum,const ActuatorVector &maximum)
{
 (void)minimum;(void)maximum;
 if(matrix_index!=0){return;}
 setSafeOutput(out);
 _status_sub.update(&_status);
 const auto now=hrt_absolute_time();
 const float dt=_last_time!=0 && now>=_last_time ? float(now-_last_time)*1e-6f : .01f;
 _last_time=now;
 const bool armed=_status.arming_state==vehicle_status_s::ARMING_STATE_ARMED;
 if(armed && !_armed_last){
  _allow_run=_run_requested && _config_valid;
  _arm_time=now;_previous={};
 }
 if(!armed){_allow_run=false;_arm_time=0;_previous={};}
 _armed_last=armed;
 // Disabling RUN stops automatic output immediately; re-enable only after disarming.
 if(!_run_requested){_allow_run=false;}
 bool finite=true;for(int i=0;i<NUM_AXES;i++){finite=finite && std::isfinite(request(i));}
 if(!armed || !_allow_run || !_config_valid || !finite
    || _status.vehicle_type!=vehicle_status_s::VEHICLE_TYPE_ROTARY_WING){
  putResidual(request,{});return;
 }
 const float elapsed=now>=_arm_time ? float(now-_arm_time)*1e-6f : 0.f;
 const float tmax=dvh::lookup(_cfg.omega_max).thrust;
 const auto scale=dvh::scales(_cfg);
 if(elapsed<_spool_seconds){
  // Keep the normal motor spool ramp, but do not leave the thrust-offset
  // pitch moment uncompensated. The same quasi-static model used by the
  // allocator predicts My = pivot_x * sum(T) + cyclic_pitch * sum(T) * u_c
  // at neutral nacelle angles. Solve that relation for the cyclic command,
  // including the current pitch-torque request.
  const float omega=_idle_omega*dvh::bound(elapsed/_spool_seconds,0.f,1.f);
  const float rotor_thrust=dvh::lookup(omega).thrust;
  const float motor_fraction=rotor_thrust/tmax;
  dvh::C5 spool_controls{{motor_fraction,motor_fraction,0.f,0.f,0.f}};
  const float total_thrust=2.f*rotor_thrust;

  if(total_thrust>1.f){
   const float requested_pitch_torque=request(1)*scale[1];
   const float pitch_bias=_cfg.pivot_x*total_thrust;
   const float pitch_per_cyclic=_cfg.cyclic_pitch*total_thrust;
   const float cyclic_ff=(requested_pitch_torque-pitch_bias)/pitch_per_cyclic;
   spool_controls[4]=dvh::bound(cyclic_ff,-1.f,1.f);
  }

  // Apply the configured cyclic slew limit so the commanded swash position
  // remains continuous when the motors leave spool-up and normal allocation
  // takes over. The physical SDF joint limits still apply independently.
  dvh::limitServoRate(_cfg,dt,_previous,spool_controls);
  _previous=spool_controls;

  out(0)=omega/_cfg.omega_max;out(1)=out(0);
  out(3+3)=_previous[4];
  out(3+6)=_previous[4];
  putResidual(request,dvh::forward(_cfg,_previous));return;
 }
 const dvh::V5 target{{request(0)*scale[0],request(1)*scale[1],request(2)*scale[2],
   request(3)*scale[3],-request(5)*scale[4]}};
 auto result=dvh::allocate(_cfg,target,_idle_omega,&_previous);
 if(!result.valid){putResidual(request,{});return;}
 dvh::limitServoRate(_cfg,dt,_previous,result.controls);
 _previous=result.controls;
 out(0)=dvh::bound(dvh::omegaForThrust(_previous[0]*tmax)/_cfg.omega_max,0.f,1.f);
 out(1)=dvh::bound(dvh::omegaForThrust(_previous[1]*tmax)/_cfg.omega_max,0.f,1.f);
 out(2)=0.f;
 out(3+3)=_previous[4];out(3+6)=_previous[4];
 out(3+5)=(_cfg.tilt_trim[0]+_previous[2]*_cfg.tilt_limit)/(dvh::pi/2.f);
 out(3+7)=(_cfg.tilt_trim[1]+_previous[3]*_cfg.tilt_limit)/(dvh::pi/2.f);
 putResidual(request,dvh::forward(_cfg,_previous));
}
void ActuatorEffectivenessDualVtol::getUnallocatedControl(int matrix_index,control_allocator_status_s &status)
{
 if(matrix_index!=0){return;}
 for(int i=0;i<3;i++){status.unallocated_torque[i]=_residual[i];status.unallocated_thrust[i]=_residual[3+i];}
}
