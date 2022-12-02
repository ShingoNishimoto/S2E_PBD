#define _USE_MATH_DEFINES
#include "RFSystemReceiver.h"
#include<cmath>

RFSystemReceiver::RFSystemReceiver(
  const int prescaler,
  ClockGenerator* clock_gen,
  const Vector<3>& compo_position_b,
  const Quaternion& q_b2c,
  Antenna* ant,
  PBD_InterSatComm* pbd_inter_sat_comm,
  const Dynamics* dynamics,
  const double update_interval_sec)
  : ComponentBase(prescaler, clock_gen),
  compo_position_b_(compo_position_b), q_b2c_(q_b2c), ant_(ant),
  dynamics_(dynamics), update_interval_sec_(update_interval_sec)
{
  fill_up(line_of_sight_c_, 0.0); line_of_sight_c_[0] = 1.0; //set(1,0,0)
  rf_system_beam_ = pbd_inter_sat_comm->GetRFSystemBeam();
}

RFSystemReceiver::~RFSystemReceiver()
{
}

void RFSystemReceiver::MainRoutine(int count)
{
  received_power_watt_ = CalcReceivedPower();
  UpdateLinkState();
}

double RFSystemReceiver::CalcReceivedPower()
{
  double power_sum = 0.0;//Sum of light intensity of each microelement[W]

  //set aparture center position
  Quaternion q_i2b = dynamics_->GetAttitude().GetQuaternion_i2b();
  Quaternion q_b2i = q_i2b.conjugate();
  Vector<3> aparture_center_pos_i = dynamics_->GetOrbit().GetSatPosition_i() + q_b2i.frame_conv(compo_position_b_);

  //line of sight in inertial coodinate system
  Quaternion q_c2b = q_b2c_.conjugate();
  Vector<3> line_of_sight_b = q_c2b.frame_conv(line_of_sight_c_);
  Vector<3> line_of_sight_i = q_b2i.frame_conv(line_of_sight_b);

  ////Find out where the beam hits the surface on which the photo detector is attached
  //Vector<3> aparture_center_from_beamwaist_i = aparture_center_pos_i - rf_system_beam_->GetBeamWaistPos_i();
  //double distance_from_beamwaist = inner_product(line_of_sight_i, aparture_center_from_beamwaist_i) / inner_product(line_of_sight_i, rf_system_beam_->GetPointingVector_i());
  //Vector<3> beam_center_on_surface_i = rf_system_beam_->GetBeamWaistPos_i() + distance_from_beamwaist * rf_system_beam_->GetPointingVector_i();

  ////If the center of the photodetector aperture and the center of the beam are significantly different, the power is considered to be zero.
  //double deviation_between_aparture_and_beam = norm(aparture_center_pos_i - beam_center_on_surface_i);
  //double beam_radius_on_surface = rf_system_beam_->CalcBeamWidthRadius(distance_from_beamwaist);
  //constexpr double cutoff_factor = 2.0;
  //if (deviation_between_aparture_and_beam > cutoff_factor * beam_radius_on_surface) { return 0.0; }

  //// Calculate received power
  //double r_from_beam_center_on_surface = norm(aparture_center_pos_i - beam_center_on_surface_i);
  //power_sum = rf_system_beam_->CalcIntensity(distance_from_beamwaist, r_from_beam_center_on_surface) * M_PI * r_aparture_m_ * r_aparture_m_;
  return power_sum;
}

void RFSystemReceiver::UpdateLinkState()
{
  if (state_ == RFLaserLinkState::Lost)
  {
    state_ = RFLaserLinkState::WaitingForLink;
    link_wait_counter_sec_ = 0.0;
  }
  else if(state_ == RFLaserLinkState::WaitingForLink)
  {
    state_ = RFLaserLinkState::Linked;
  }
  else //state_ == RFLaserLinkState::Linked
  {
    state_ = RFLaserLinkState::Lost;
  }
}


std::string RFSystemReceiver::GetLogHeader() const
{
  std::string str_tmp = "";
  str_tmp += WriteScalar("RF received power", "W");
  return str_tmp;
}

std::string RFSystemReceiver::GetLogValue() const
{
  std::string str_tmp = "";
  str_tmp += WriteScalar(received_power_watt_);
  return str_tmp;
}
