#include "InitGNSSReceiver.hpp"
#include "../../Simulation/PBD/InitGNSSAntennaPCC.hpp"

#include <string.h>
#include <Dense>
#include "Interface/InitInput/IniAccess.h"

typedef struct _gnssrecever_param {
  int prescaler;
  AntennaModel antenna_model;
  Vector<3> antenna_pos_b;
  Quaternion q_b2c;
  double half_width;
  std::string gnss_id;
  int ch_max;
  Vector<3> noise_std;
  Vector<3> alignment_err_std;
  std::string antex_file_path;
  double pcv_d_azimuth;
  double pcv_d_elevation;
  double pseudo_stddev;
  double carrier_stddev;
  double clock_rn_stddev; // for random noise
} GNSSReceiverParam;


GNSSReceiverParam ReadGNSSReceiverIni(const std::string fname, const GnssSatellites* gnss_satellites) {
  GNSSReceiverParam gnssreceiver_param;

  IniAccess gnssr_conf(fname);
  char GSection[30] = "GNSSReceiver";

  int prescaler = gnssr_conf.ReadInt(GSection, "prescaler");
  if (prescaler <= 1) prescaler = 1;
  gnssreceiver_param.prescaler = prescaler;

  gnssreceiver_param.antenna_model = static_cast<AntennaModel>(gnssr_conf.ReadInt(GSection, "antenna_model"));
  if (!gnss_satellites->IsCalcEnabled() && gnssreceiver_param.antenna_model == CONE) {
    std::cout << "Calculation of GNSS SATELLITES is DISABLED, so the antenna "
                 "model of GNSS Receiver is automatically set to SIMPLE model."
              << std::endl;
    gnssreceiver_param.antenna_model = SIMPLE;
  }

  gnssr_conf.ReadVector(GSection, "antenna_pos_b", gnssreceiver_param.antenna_pos_b);
  gnssr_conf.ReadQuaternion(GSection, "q_b2c", gnssreceiver_param.q_b2c);
  gnssreceiver_param.half_width = gnssr_conf.ReadDouble(GSection, "half_width");
  gnssreceiver_param.gnss_id = gnssr_conf.ReadString(GSection, "gnss_id");
  gnssreceiver_param.ch_max = gnssr_conf.ReadInt(GSection, "ch_max");
  gnssr_conf.ReadVector(GSection, "nr_stddev_eci", gnssreceiver_param.noise_std);
  gnssr_conf.ReadVector(GSection, "alignment_err_stddev_b", gnssreceiver_param.alignment_err_std);
  const std::string antex_file_name = gnssr_conf.ReadString(GSection, "antex_file_name");
  // ここには必要ないかもしれない．FIXME
  if (antex_file_name.substr(antex_file_name.length() - 3, 3) == "atx")
  {
    gnssreceiver_param.antex_file_path = "../../../ExtLibraries/ANTEX/" + antex_file_name;
  }
  else if (antex_file_name.substr(antex_file_name.length() - 3, 3) == "csv")
  {
    gnssreceiver_param.antex_file_path = antex_file_name; // こっちはパスを直接指定しておく．
  }
  else
  {
    std::cout << "AntexFileError!" << std::endl;
  }
  gnssreceiver_param.pcv_d_azimuth = gnssr_conf.ReadDouble(GSection, "d_azi");
  gnssreceiver_param.pcv_d_elevation = gnssr_conf.ReadDouble(GSection, "d_ele");
  gnssreceiver_param.pseudo_stddev = gnssr_conf.ReadDouble(GSection, "pseudo_range_stddev");
  gnssreceiver_param.carrier_stddev = gnssr_conf.ReadDouble(GSection, "carrier_phase_stddev") / 1000.0; // mm -> m
  gnssreceiver_param.clock_rn_stddev = gnssr_conf.ReadDouble(GSection, "clock_rn_stddev");

  return gnssreceiver_param;
}

PBD_GNSSReceiver InitGNSSReceiver(ClockGenerator* clock_gen, int id, const std::string fname, const Dynamics* dynamics,
                              const GnssSatellites* gnss_satellites, const SimTime* simtime)
{
  GNSSReceiverParam gr_param = ReadGNSSReceiverIni(fname, gnss_satellites);

  PhaseCenterCorrection* pcc = InitPCC(gr_param.antex_file_path, gr_param.pcv_d_azimuth, gr_param.pcv_d_elevation);


  PBD_GNSSReceiver gnss_r(gr_param.prescaler, clock_gen, id, gr_param.gnss_id, gr_param.ch_max, gr_param.antenna_model, gr_param.antenna_pos_b,
                      gr_param.q_b2c, gr_param.half_width, gr_param.noise_std,
                      gr_param.alignment_err_std, pcc,
                      gr_param.pcv_d_azimuth, gr_param.pcv_d_elevation,
                      gr_param.pseudo_stddev, gr_param.carrier_stddev, gr_param.clock_rn_stddev,
                      dynamics, gnss_satellites, simtime);
  return gnss_r;
}
