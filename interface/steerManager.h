#ifndef __STEERMAN__
#define __STEERMAN__

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <list>

class steerManager {
   public :
      steerManager();
      steerManager(std::string _fname_);
      ~steerManager();

      void setFname(std::string _fname_);
      void addSingleParameter(std::string _singleParName_);
      void addArrayParameter(std::string _arrayParName_);

      bool read();

      void printPars();

      float *getArrayPara(std::string _parname_ );
      float getSinglePara(std::string _parname_ );

      float getCorrectionAtPoint(float _value_, std::string _edgesName_,
            std::string _corrArrayName_);

   private :
      std::string _steerFileName;

      std::list < std::string > singleParamList;
      std::list < std::string > arrayParamList;

      std::map < std::string,std::list < float > > _vectorParas;
      std::map < std::string,float > _singleParas;
};
#endif
