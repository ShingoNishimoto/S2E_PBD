#include "RFSystemTransmitter.h"
#include "Matrix.hpp"
#include "MatVec.hpp"
#include "NormalRand.hpp"
#include "GlobalRand.h"

RFSystemTransmitter::RFSystemTransmitter(
  int prescaler,
  ClockGenerator* clock_gen,
  const Vector<3>& compo_position_b,
  const Quaternion& q_b2c,
  ANT* ant,
  PBD_InterSatComm* pbd_inter_sat_comm,
  const Dynamics* dynamics
) : ComponentBase(prescaler, clock_gen),
compo_position_b_(compo_position_b), q_b2c_(q_b2c), ant_(ant), 
dynamics_(dynamics)
{
  fill_up(line_of_sight_c_, 0.0); line_of_sight_c_[0] = 1.0; //(1,0,0)
  rf_system_beam_ = pbd_inter_sat_comm->GetRFSystemBeam();
}

RFSystemTransmitter::~RFSystemTransmitter()
{
}

void RFSystemTransmitter::MainRoutine(int count)
{
  //calculate raw diffraction angle
  //emission_angle_rad_[0] = (rf_system_beam_->GetWaveLength() * acoustic_signal_freq_hz_[0]) / (2.0 * acoustic_signal_speed_in_lattice_m_s_);
  //emission_angle_rad_[1] = (rf_system_beam_->GetWaveLength() * acoustic_signal_freq_hz_[1]) / (2.0 * acoustic_signal_speed_in_lattice_m_s_);

  ////set frequency shift
  //rf_system_beam_->SetFreqShift(2.0 * acoustic_signal_freq_hz_);

  //set position of beam waist
  Quaternion q_i2b = dynamics_->GetAttitude().GetQuaternion_i2b();
  Quaternion q_b2i = q_i2b.conjugate();
  Vector<3> pos_beamwaist_i = dynamics_->GetOrbit().GetSatPosition_i() + q_b2i.frame_conv(compo_position_b_);
  //rf_system_beam_->SetBeamWaistPos_i(pos_beamwaist_i);

  AddNoise();

  ////set beam pointing
  //Vector<3> pointing_vec_c = CalcPointingVector_c();
  /*Quaternion q_c2b = q_b2c_.conjugate();
  Vector<3> pointing_vec_b = q_c2b.frame_conv(pointing_vec_c);
  Vector<3> pointing_vec_i = q_b2i.frame_conv(pointing_vec_b);*/
  //rf_system_beam_->SetPointingVector_i(pointing_vec_i); //inertial coordinate system

  // Adjust optical params with beam expander
  ExpandBeam();
}

void RFSystemTransmitter::ExpandBeam()
{
  // [TODO] Implement
}

//Vector<3> RFSystemTransmitter::CalcPointingVector_c()
//{
//  Matrix<3, 3> rotation_mat_theta1;
//  rotation_mat_theta1[0][0] = cos(emission_angle_rad_[0]); rotation_mat_theta1[0][1] = -sin(emission_angle_rad_[0]); rotation_mat_theta1[0][2] = 0.0;
//  rotation_mat_theta1[1][0] = sin(emission_angle_rad_[0]); rotation_mat_theta1[1][1] =  cos(emission_angle_rad_[0]); rotation_mat_theta1[1][2] = 0.0;
//  rotation_mat_theta1[2][0] = 0.0;                         rotation_mat_theta1[2][1] =  0.0;                         rotation_mat_theta1[2][2] = 1.0;
//
//  Matrix<3, 3> rotation_mat_theta2;
//  rotation_mat_theta2[0][0] = cos(emission_angle_rad_[1]);  rotation_mat_theta2[0][1] = 0.0;                         rotation_mat_theta2[0][2] = sin(emission_angle_rad_[1]);
//  rotation_mat_theta2[1][0] = 0.0;                          rotation_mat_theta2[1][1] = 1.0;                         rotation_mat_theta2[1][2] = 0.0;
//  rotation_mat_theta2[2][0] = -sin(emission_angle_rad_[1]); rotation_mat_theta2[2][1] = 0.0;                         rotation_mat_theta2[2][2] = cos(emission_angle_rad_[1]);
//
//  Matrix<3, 3> rotation_matrix = rotation_mat_theta1 * rotation_mat_theta2;
//  Vector<3> pointing_vec_c = rotation_matrix * line_of_sight_c_;
//
//  return pointing_vec_c;
//}

void RFSystemTransmitter::AddNoise()
{
  // add pointing noise
  
  // add emission power noise
  //double raw_total_power = rf_system_beam_->GetTotalPower();
  //double emission_power_stddev = raw_total_power * emission_power_stddev_ratio_;
  //libra::NormalRand power_with_noise(raw_total_power, emission_power_stddev, g_rand.MakeSeed());

  //rf_system_beam_->SetTotalPower((double)power_with_noise);
}

string RFSystemTransmitter::GetLogHeader() const
{
  string str_tmp = "";
  //str_tmp += WriteVector("acoustic_signal_freq_center", "", "Hz", 2);
  return str_tmp;
}

string RFSystemTransmitter::GetLogValue() const
{
  string str_tmp = "";
  //str_tmp += WriteVector(acoustic_signal_freq_center_hz_, 9);
  return str_tmp;
}
