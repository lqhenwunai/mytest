#ifndef __Message__
#define __Message__

//#include <iostream>
//#include <fstream>
//#include <iomanip>
//#include <cstdio>
//#include <stdlib.h>
//

class Message
{
public:
   void Header(const char* msg);
   void ConvFail(const char* msg);
   void FileError(const char* msg);
   void PointerError(const char* msg);
   void DimError(const char* header);

};


#endif
