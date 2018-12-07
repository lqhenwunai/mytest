#include "Message.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <stdlib.h>
using namespace std;
void Message::Header(const char* header)
{
   cout<<"-----------------------------------------------"<<endl;
   cout<<"             "<<header<<"                   "<<endl;
   cout<<"-----------------------------------------------"<<endl;

}


void Message::ConvFail(const char* header)
{
   cout<<"-----------------------------------------------"<<endl;
   cout<<"  Convergence failed:         "<<header<<"                   "<<endl;
   cout<<"-----------------------------------------------"<<endl;
   abort();
}

void Message::DimError(const char* header)
{
   cout<<"-----------------------------------------------"<<endl;
   cout<<"  Dimension error      "<<header<<"            "<<endl;
   cout<<"-----------------------------------------------"<<endl;
   abort();
}

void Message::FileError(const char* header)
{
   cout<<"-----------------------------------------------"<<endl;
   cout<<"  file open error:         "<<header<<"                   "<<endl;
   cout<<"-----------------------------------------------"<<endl;
   abort();
}

