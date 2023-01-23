#include "InitGNSSAntennaPCC.hpp"

#include <string.h>
#include <Dense>
#include "Interface/InitInput/IniAccess.h"

using std::vector;
using std::string;

static bool ReadAntexTable(string file_name, const double d_azi, const double d_ele, libra::Vector<3>& pco, std::vector<double>& pcv);
static bool ReadAntexFromCsv(string file_name, const double d_azi, const double d_ele, libra::Vector<3>& pco, std::vector<double>& pcv);
static vector<string> split(string& input, char delimiter);

PhaseCenterCorrection* InitPCC(string file_name, const double d_azi, const double d_ele)
{
  libra::Vector<3> pco_;
  vector<double> pcv_;

  if (file_name.substr(file_name.length() - 3, 3) == "atx")
  {
    if (!ReadAntexTable(file_name, d_azi, d_ele, pco_, pcv_))
    {
      // モデルを球面PCO 0，PCV球面のモデルにする？
      abort();
    }
  }
  else if (file_name.substr(file_name.length() - 3, 3) == "csv")
  {
    if (!ReadAntexFromCsv(file_name, d_azi, d_ele, pco_, pcv_))
    {
      abort();
    }
  }
  else
  {
    std::cout << "AntexFileError!" << std::endl;
  }

  PhaseCenterCorrection* pcc = new PhaseCenterCorrection(pco_, pcv_, d_azi, d_ele);
  return pcc;
}


static bool ReadAntexTable(string file_name, const double d_azi, const double d_ele, libra::Vector<3>& pco, std::vector<double>& pcv)
{
  std::ifstream antex_file(file_name);
  if (!antex_file.is_open()) {
    std::cerr << "file open error:Antex\n";
    return false;
  }

  string line;
  // PCO
  std::getline(antex_file, line);
  std::istringstream streamline(line);
  // string description;
  streamline >> pco[1] >> pco[0] >> pco[2]; // north, east, upなのでy, x, zの順で入れる．
  // streamline >> pco[0] >> pco[1] >> pco[2];
  // pco[1] *= -1; // x, y, zの順で代入した場合はyはwestなので符号を反転．

  // NOAZI (skip)
  std::getline(antex_file, line);

  // PCV
  int num_azimuth = (360 / d_azi + 1);
  int num_elevation = (90 / d_ele + 1);
  int num_values = num_azimuth * num_elevation;
  for (int i = 0; i < num_azimuth; i++) {
    string line;
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

static bool ReadAntexFromCsv(string file_name, const double d_azi, const double d_ele, libra::Vector<3>& pco, std::vector<double>& pcv)
{
  std::ifstream csv_file(file_name);
  if (!csv_file.is_open()) {
    std::cerr << "file open error:Csv\n";
    return false;
  }

  string line;
  // PCO
  std::getline(csv_file, line);
  vector<string> str_pco_vec = split(line, ',');
  for (int i = 0; i < 3; i++) pco[i] = std::stod(str_pco_vec[i]); // logに入っているのものなのでx, y, zの順で入れる．

  // PCV
  int num_azimuth = (360 / d_azi + 1);
  int num_elevation = (90 / d_ele + 1);
  int num_values = num_azimuth * num_elevation;
  for (int i = 0; i < num_azimuth; i++) {
    string line;
    std::getline(csv_file, line);
    vector<string> str_pcv_vec = split(line, ',');

    for (int j = 0; j < num_elevation; j++)
    {
      double pcv_element = std::stod(str_pcv_vec[j]);
      pcv.push_back(pcv_element);
    }
  }
  return true;
}

static vector<string> split(string& input, char delimiter)
{
    std::istringstream stream(input);
    string field;
    vector<string> result;
    while (std::getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}
