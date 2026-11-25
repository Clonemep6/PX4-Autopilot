#pragma once
// SITL research prototype. See package README before enabling motor output.
#include "control_allocation/actuator_effectiveness/ActuatorEffectiveness.hpp"
#include <px4_platform_common/module_params.h>
#include <uORB/Subscription.hpp>
#include <uORB/topics/vehicle_status.h>
#include <drivers/drv_hrt.h>
#include <lib/dual_vtol_hover/DualVtolHoverMath.hpp>

class ActuatorEffectivenessDualVtol : public ModuleParams, public ActuatorEffectiveness
{
public:
 explicit ActuatorEffectivenessDualVtol(ModuleParams *parent);
 bool getEffectivenessMatrix(Configuration &configuration, EffectivenessUpdateReason reason) override;
 void updateSetpoint(const matrix::Vector<float, NUM_AXES> &control_sp, int matrix_index,
                     ActuatorVector &actuator_sp, const ActuatorVector &actuator_min,
                     const ActuatorVector &actuator_max) override;
 void getUnallocatedControl(int matrix_index, control_allocator_status_s &status) override;
 const char *name() const override { return "Dual VTOL hover (experimental SITL)"; }
private:
 void updateParams() override;
 void putResidual(const matrix::Vector<float, NUM_AXES> &request, const dvh::V5 &achieved);
 void setSafeOutput(ActuatorVector &out) const;
 float readFloat(const char *name, float fallback);
 int32_t readInt(const char *name, int32_t fallback);
 dvh::Config _cfg{};
 dvh::C5 _previous{};
 float _idle_omega{400.f};
 float _spool_seconds{5.f};
 bool _run_requested{false};
 bool _dvt_enabled{false};
 float _dvt_push_test{0.f};
 bool _dvt_ramp_enabled{false};
 float _dvt_pusher_output{0.f};
 bool _config_valid{false};
 bool _armed_last{false};
 bool _allow_run{false};
 bool _reported_invalid{false};
 hrt_abstime _arm_time{0};
 hrt_abstime _last_time{0};
 float _residual[NUM_AXES]{};
 uORB::Subscription _status_sub{ORB_ID(vehicle_status)};
 vehicle_status_s _status{};
};
