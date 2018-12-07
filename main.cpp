//#ifndef __main__
//#define __main__
 
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <stdlib.h>
#include "test.h" 
#include "FileRead.h"
#include "Matrix-tools.h"
#include "Message.h"
#include "Util.h"
#include "Matrix-math.h"
//#endif
//----------------------
//on input:  C, MO coef.
//           occ, 
//on output: D
//--------------------

int Density_calc(TRMatrix &D, TRMatrix &C, int occ)
{
  double tmp;
  for(int i=0; i<C.rows; i++){
    for(int j=0; j<C.cols; j++){
       tmp=0.0;
       for(int m=0; m<occ; m++){  //10 electrons for water, 5 orbitals
         tmp+=C[i][m]*C[j][m];
      }
       D[i][j]=tmp;
    }
  }
 return 0;
}

//----------------------
//Diagnoliize fock matrix
//on input:  S12, S^(-1/2)
//           Fold, old Fock matrix 
//on output: Fock
//           C0, MO coef
//--------------------
 
int Fock_calc(TRMatrix &Fock, TRMatrix &S12, TRMatrix& Fold, TRMatrix &C0)
{
  int row,col,res;
  row=Fold.rows;
  col=Fold.cols;
  TRMatrix ST;
  ST.NewMat(row,col);
  ST.Transpose(S12);
  Mat_x_Mat_x_Mat(Fock,ST,Fold,S12);//actually ST==inv_sqrt_overlap

//------------------------------
//diagonialize the Fock matrix
//obtain the eigenvector in orthonormal AO basis
//-------------------------------  
  TRMatrix Cprime;
  Cprime.NewMat(row,col);
  Cprime.Unity();  
  res=Fock.Diagnolization(Cprime,1e-14);//each coloumn is an eigenvector corresponding to the eigenvalue, and mutually orthogonal to each other

  Math_Mat_x_Mat(C0,S12,Cprime);//C=S^(-1/2)*C'
  return 0;
}



using namespace std;
int main()
{
//   int iter=125;//scf iteraction

   Message msg;
   ifstream infile;
   infile.open("overlap.dat");
   if(!infile)
   {  msg.FileError("overlap/dat");}

   int res,count=0,row=0, col=0;
   int* mu, *nu;
   double *value;
   double threshold_diag=1e-14, tmp=0.0;

   count=GetLineCount(infile);  
   mu = new int[count];
   nu = new int[count];
   value = new double[count];
//---------------
//overlap matrix
//--------------
   TRMatrix overlap; 
   FileRead(infile, count, row, col,mu,nu,value);
   overlap.NewMat(row,col);  //read the low triangle and make it full matrix
   overlap.ReadMat(value);
   overlap.msg.Header("overlap matrix");
   cout<<"overlap matrix"<<endl;
   print_matrix(overlap.p, row, col);
   cout<<endl;

   TRMatrix overlap_inv;
   overlap_inv.NewMat(row, col);
   overlap_inv.CopyMat(overlap);
   res=Inverse_Mat(overlap_inv);

//   TRMatrix Id;
//   Id.NewMat(row,col);
//   Math_Mat_x_Mat(Id, overlap, overlap_inv);
//    
//   cout<<"Id matrix"<<endl;
//   print_matrix(Id);


//    TRMatrix U,L;
//    U.NewMat(row, col);
//    L.NewMat(row, col);
//    LU_Decomposition(overlap,L,U);
//    cout<<"L"<<endl;
//    print_matrix(L);
//    cout<<"U"<<endl;
//    print_matrix(U);
//--------------------------
//obtain S^(-1/2)
//-------------------------

   TRMatrix EigenVector, EigenVector_T;
   EigenVector.NewMat(row,col);
   EigenVector_T.NewMat(row,col);
   EigenVector.Unity();
   cout<<"start diagnolization"<<endl;

//obtain diagnoalized overlap_inv  S^-1   
   res=overlap_inv.Diagnolization(EigenVector,threshold_diag);
   if(res==1)
    {overlap.msg.Header("overlap matrix diagnolized");
     print_matrix(overlap_inv);
    }
   else
    {overlap.msg.ConvFail("overlap diagnolization");}
   
//get the square root of diagnolized overlp_inv
   for(int i=0; i<overlap_inv.rows; i++){
      overlap_inv[i][i]=sqrt(overlap_inv[i][i]);
   }
   cout<<"diagnolized overlap inverse square root"<<endl; 
   print_matrix(overlap_inv);
  
   EigenVector_T.Transpose(EigenVector);

   TRMatrix inv_sqrt_overlap;
   inv_sqrt_overlap.NewMat(row,col);

//   cout<<"eigen vectors"<<endl;
//   print_matrix(EigenVector);
//   cout<<"eigen vectors Transpose"<<endl;
//   print_matrix(EigenVector_T);

// S^(-1/2)!=its transpose
//http://www.psicode.org/psi4manual/1.2/scf.html#orthogonalization
//http://booksite.elsevier.com/9780444594365/downloads/16755_10030.pdf
   cout<<"symmetric orthogonalization matrix S^(-1/2)"<<endl;
   Mat_x_Mat_x_Mat(inv_sqrt_overlap, EigenVector,overlap_inv,EigenVector_T);
   print_matrix(inv_sqrt_overlap);


   EigenVector.DelMat();
   EigenVector_T.DelMat();
   overlap_inv.DelMat();
//-----------------
//kinetic matrix
//----------------
   
   TRMatrix kinetic;
   infile.open("kinetic.dat"); 
   FileRead(infile, count, row, col,mu,nu,value);
   kinetic.NewMat(row,col);
   kinetic.ReadMat(value);


//------------------
//nuclear attraction
//------------------
   TRMatrix nuc_att;
   infile.open("Nuc.dat");
   FileRead(infile, count, row, col,mu,nu,value);
   nuc_att.NewMat(row,col);
   nuc_att.ReadMat(value);


   TRMatrix Hcore;
   Hcore.NewMat(row, col);
   for (int i=0; i<row; i++){
     for (int j=0; j<col; j++){
         Hcore[i][j]=kinetic[i][j]+nuc_att[i][j];
     }
    }

   msg.Header("Hcore matrix");
   Hcore.WriteMat();

   delete []mu,nu,value;

//-------------------------------------------------
//form guess fock matrix in the orthonomal AO basis
//F0=ST^(-1/2)*Hcore*S^(-1/2)
//FC=eSC
//(S-1/2)*F*(S-1/2)*(S1/2)*C=e(S-1/2)*S*C=eS(1/2)*C
//------------------------------------------------
  TRMatrix F0,ST,C0;//Fock matrix, transpose of S(-1/2), MO coeff
  F0.NewMat(row,col);
  ST.NewMat(row,col);
  C0.NewMat(row,col);
  Fock_calc(F0,inv_sqrt_overlap,Hcore,C0);

  msg.Header("Initial Fock Matirx in the orthonormal AO basis");
  print_matrix(F0);

  msg.Header("diagnolized F0 Matrix");
  print_matrix(F0);
  cout<<"Initial MO coefficient"<<endl;
  print_matrix(C0);

//----------------------
//density
//------------------------
  msg.Header("initial density");
  TRMatrix D;//density
  D.NewMat(row,col);

  int occ=5;//10 electrons for water, 5 orbitals
  Density_calc(D,C0,occ);
  print_matrix(D);

//------------------
//initial SCF energy
//------------------
  double Eele=0.0;
  for(int mu=0; mu<row; mu++){
    for(int nu=0; nu<col; nu++){
       Eele+=D[mu][nu]*(Hcore[mu][nu]*2.0);
    }
  }
    
  double Etotal=0.0;
  double ENN=8.002367061810450; 
  Etotal=Eele+ENN;
  msg.Header("SCF iterations");
  cout<<"Iter          "<<"E(elec)          "<<"Etot"<<endl;
  cout<<"00  "<<Eele<<"     "<<Etotal<<endl;
//  print_matrix(F0);
//  cout<<Eele2<<endl;
  double Eold=Eele;  
//--------------------------
//two electron integrals
//--------------------------
  infile.open("two-electron.dat");  
  count=GetLineCount(infile);

  int  *rho, *lemda;
  mu    = new int[count];
  nu    = new int[count];
  rho   = new int[count];
  lemda = new int[count];
  value = new double[count];
  int MU, NU, RHO, LEMDA; //dimesion of pointers
  
  FileRead_2(infile, count, MU, NU, RHO, LEMDA,mu,nu,rho,lemda,value);
  TRMatrixContainer I;
  I.NewMat(MU, NU, RHO, LEMDA);
  I.ReadMat(value,mu,nu,rho,lemda, count);

//  print_matrixcontainer(I);
  

//---------------------------
//compute new fock matrix
//--------------------------
  TRMatrix Fock;
  Fock.NewMat(row,col);
  print_matrix(Fock);
  print_matrix(D);
  for(int mu=0; mu<row; mu++){
   for(int nu=0; nu<col; nu++){
     Fock[mu][nu]=Hcore[mu][nu];
     for(int lemda=0; lemda<row; lemda++){
       for(int sigma=0; sigma<col; sigma++){
          Fock[mu][nu]+=D[lemda][sigma]*(2.0*I(mu,nu,lemda,sigma)-I(mu,lemda,nu,sigma));
       }//sigma
      }//lemda
    }//nu
   }//mu

  msg.Header("new fock matrix");
  print_matrix(Fock);

//--------------------
//new density matrix
//Fock'=ST^(-1/2)Fock*S^(-1/2)
//--------------------

  TRMatrix Fnew,Cnew,Dnew,Fnew_unDiag;
  Fnew.NewMat(row,col);
  Cnew.NewMat(row,col);
  Dnew.NewMat(row,col);
  Fnew_unDiag.NewMat(row,col);//not diagonalized fock matrix.
  Fnew_unDiag.CopyMat(Fock);
  Fock_calc(Fnew,inv_sqrt_overlap,Fock,Cnew);
//  cout<<"cnew"<<endl;
//  print_matrix(Cnew);
//  print_matrix(inv_sqrt_overlap);
  Density_calc(Dnew,Cnew,occ);

//  print_matrix(Fnew); 
//---------------------
//new scf energy
//---------------------
 Eele=0.0;
  for(int mu=0; mu<row; mu++){
    for(int nu=0; nu<col; nu++){
       Eele+=Dnew[mu][nu]*(Hcore[mu][nu]+Fnew_unDiag[mu][nu]);
    }
  }

  Etotal=Eele+ENN;

  cout<<"Iter          "<<"E(elec)          "<<"Etot"<<endl;
  cout<<"01  "<<Eele<<"     "<<Etotal<<endl;


  double Enew=Eele;
  double deltaE=1e-8;
  double RMSD=1e-8;
  double tmpE, tmpD=0.0;
  int iter=0;
/*
  while (iter<125)
{
  iter++;
  tmpE=Eele-Eold;
  for(int mu=0; mu<occ; mu++){
    for(int nu=0; nu<occ; nu++){
       tmpD+=(Dnew[mu][nu]-D[mu][nu])*(Dnew[mu][nu]-D[mu][nu]);
    }
  }
  tmpD=sqrt(tmpD);

  if(tmpE<deltaE and tmpD<RMSD)
    { 
      msg.Header("SCF converged");
      break;
    }
//construct new Fock matrix
  for(int mu=0; mu<row; mu++){
   for(int nu=0; nu<col; nu++){
     Fock[mu][nu]=Hcore[mu][nu];
     for(int lemda=0; lemda<row; lemda++){
       for(int sigma=0; sigma<col; sigma++){
          Fock[mu][nu]+=Dnew[lemda][sigma]*(2.0*I(mu,nu,lemda,sigma)-I(mu,lemda,nu,sigma));
       }//sigma
      }//lemda
    }//nu
   }//mu


//--------------------
//diagnolize Fock
//new density matrix
//--------------------

  Cnew.Unity();
  Dnew.Init();

  Fock_calc(Fnew,inv_sqrt_overlap,Fock,Cnew);
  Density_calc(Dnew,Cnew,occ);

//---------------------
//new scf energy
//---------------------
 Eele=0.0;
  for(int mu=0; mu<row; mu++){
    for(int nu=0; nu<col; nu++){
       Eele+=Dnew[mu][nu]*(Hcore[mu][nu]+Fock[mu][nu]);
    }
  }

  Etotal=Eele+ENN;


}
  

*/
 
  return 0;
}


