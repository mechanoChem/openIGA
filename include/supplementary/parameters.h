/**
* Provide parameters storage
*/

#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#include <iostream>
#include <sstream>
#include <string.h>
#include <map>
#include <deal.II/base/point.h>
class parametersClass{
 public:
	 /**
	 *
	 **/
  void setDouble(std::string param, double value, bool print=false);
  void setInt(std::string param, int value, bool print=false);
  void setBool(std::string param, bool value, bool print=false);
  void setString(std::string param, std::string value, bool print=false);
	
	void setPoint(std::string param, dealii::Point<3> value, bool print=false);
	void setPoint(std::string param, dealii::Point<2> value, bool print=false);
	void setPoint(std::string param, dealii::Point<1> value, bool print=false);
	
  double getDouble(std::string param);
  int getInt(std::string param);
  bool getBool(std::string param);
  std::string getString(std::string param);
	dealii::Point<3> getPoint(std::string param);
	dealii::Point<2> getPoint(std::string param);
	dealii::Point<1> getPoint(std::string param);
  //void readInParameters(std::string fileName);
 private:
  std::map<std::string, int> pInt;
  std::map<std::string, double> pDouble;
  std::map<std::string, bool> pBool;
  std::map<std::string, std::string> pString;
	std::map<std::string, dealii::Point<3> > pPoint3;
	std::map<std::string, dealii::Point<2> > pPoint2;
	std::map<std::string, dealii::Point<1> > pPoint1;
};

#endif