#include "InitGNSSReceiver.hpp"

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
} GNSSReceiverParam;

bool ReadAntexTable(std::string file_name, const double d_azi, const double d_ele, libra::Vector<3>& pco, std::vector<double>& pcv);

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
  gnssreceiver_param.antex_file_path = "../../../ExtLibraries/ANTEX/" + gnssr_conf.ReadString(GSection, "antex_file_name");
  gnssreceiver_param.pcv_d_azimuth = gnssr_conf.ReadDouble(GSection, "d_azi");
  gnssreceiver_param.pcv_d_elevation = gnssr_conf.ReadDouble(GSection, "d_ele");

  return gnssreceiver_param;
}

PBD_GNSSReceiver InitGNSSReceiver(ClockGenerator* clock_gen, int id, const std::string fname, const Dynamics* dynamics,
                              const GnssSatellites* gnss_satellites, const SimTime* simtime)
{
  GNSSReceiverParam gr_param = ReadGNSSReceiverIni(fname, gnss_satellites);

  libra::Vector<3> pco_(0);
  std::vector<double> pcv_;

  if (!ReadAntexTable(gr_param.antex_file_path, gr_param.pcv_d_azimuth, gr_param.pcv_d_elevation, pco_, pcv_))
  {
    // hoge???
    // モデルを球面PCO 0，PCV球面のモデルにする？
  }

  PBD_GNSSReceiver gnss_r(gr_param.prescaler, clock_gen, id, gr_param.gnss_id, gr_param.ch_max, gr_param.antenna_model, gr_param.antenna_pos_b,
                      gr_param.q_b2c, gr_param.half_width, gr_param.noise_std, gr_param.alignment_err_std, pco_, pcv_, gr_param.pcv_d_azimuth, gr_param.pcv_d_elevation, dynamics, gnss_satellites, simtime);
  return gnss_r;
}

bool ReadAntexTable(std::string file_name, const double d_azi, const double d_ele, libra::Vector<3>& pco, std::vector<double>& pcv)
{
  std::ifstream antex_file(file_name);
  if (!antex_file.is_open()) {
    std::cerr << "file open error:Antex\n";
    return false;
  }

  std::string line;
  // PCO
  std::getline(antex_file, line);
  std::istringstream streamline(line);
  // std::string description;
  streamline >> pco[0] >> pco[1] >> pco[2];

  // NOAZI (skip)
  std::getline(antex_file, line);

  // PCV
  int num_azimuth = (360 / d_azi + 1);
  int num_elevation = (90 / d_ele + 1);
  int num_values = num_azimuth * num_elevation;
  for (int i = 0; i < num_azimuth; i++) {
    std::string line;
    std::getline(antex_file, line);
    std::istringstream streamline(line);

    double azimuth;
    streamline >> azimuth; // azimuth angle
    for (int j = 0; j < num_elevation; j++)
    {
      double pcv_element;
      streamline >> pcv_element;
      pcv.push_back(pcv_element);
    }
  }
  return true;
}
