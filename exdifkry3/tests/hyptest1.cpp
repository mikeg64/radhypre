// Copyright 2005, Google Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// A sample program demonstrating using Google C++ testing framework.

#include "margtest1.h"
#include "../../include/octant.h"
//#include "radtmesh.h"
#include "radraytsimulation.h"
#include "RadInputData.h"
#include <iostream>



//only a single octant!
double ParamFileRead()
{
	double result;
  RadInputData data;
	
  data.display();
   std::cout << "test read input file" << std::endl;
   if (data.readInputFile("input.txt")) {
        data.display();  // Display parsed data
    } else {
        std::cerr << "Failed to read input file!" << std::endl;
    }


	result=0;
	return result;
}

//The model is a single mesh composed of
// a collection of octants
// each octant will have a list of grids
// each grid will have a single field


//only a single octant!
double InstantiateOctant()
{
	double result;
	octant myoctant(3,2.0,2.0,2.0);  //create an octant with 3 grids at 2,2,2


	result=myoctant.disp[0];
	return result;
}

//call the mesh constructor 
//which will have a single octant by default
int InstantiateRadtmesh()
{
	int result;
  RadInputData data;
  radtmesh simmesh(data);  //10 timesteps, 1 level and 1 grid per level

	std::shared_ptr<grid>thisgrid; 
  std::shared_ptr<octant>thisoctant; 

  thisoctant=std::static_pointer_cast<octant>(simmesh.getoctant(0));
  if(thisoctant)
  {
     double x=thisoctant->getposx();
     cout << "\n \ninst radtmesh " << x <<  "    " << simmesh.octants.size() << '\n';
  }
  else
    cout << "\n \noctant is null "  << '\n';

	result=simmesh.m_na;
	return result;
}


//checking we can create a field and push it to a list and recover it
int CreateField()
{
  int result;
  RadInputData data;
  radtmesh simmesh(data);

  std::shared_ptr<rad>rf;
  std::shared_ptr<rad>crad =std::make_shared<rad>(simmesh.m_nf,simmesh.m_na,1);  // create shared pointers to radiation fields
  std::shared_ptr<rad>crad1 =std::make_shared<rad>(simmesh.m_nf,simmesh.m_na,1);
  std::vector<std::shared_ptr<rad >> rads;
  rf=crad1;
  crad->radtemp=30;
  //rf=std::static_pointer_cast<rad>(crad->copy());
  rf->radtemp=50;
  rads.push_back(rf);
  rads.push_back(crad);

  std::cout << "radtemp "<<crad->radtemp << std::endl;
  std::cout << "radtemp "<<rf->radtemp << std::endl;
  //rf=(crad->copy());

	result=simmesh.m_na;
	return result;  
}


int InitialiseGrid()
{
	int result;
  int i1,i2,i3;

  RadInputData data;
  radtmesh simmesh(data);


   int n1=simmesh.nx;//tg->n1;
    int n2=simmesh.ny;//tg->n2;
    int n3=simmesh.nz;//tg->n3; //number of cells in each direction

i1=0;
i2=0;
i3=0;
 std::shared_ptr<rad>rf;
	std::shared_ptr<grid>thisgrid; 
  std::shared_ptr<octant>thisoctant; 
std::shared_ptr<grid> tg=simmesh.octants[0]->mgrid[1];
  thisoctant=std::static_pointer_cast<octant>(simmesh.getoctant(0));
  if(thisoctant)
  {
     double x=thisoctant->getposx();
     cout << "\n \ninst radtmesh " << x <<  "    " << simmesh.octants.size() << '\n';
  }
  else
    cout << "\n \noctant is null "  << '\n';

  simmesh.creategrid();
  simmesh.initgrid();


  for(i1=0;i1<simmesh.m_na;++i1)
        std::cout << simmesh.m_adc[i1] <<" ";
  std::cout<<std::endl;


    //rf= std::dynamic_pointer_cast<rad>((tg->ifield[2]));
    //rf= std::static_pointer_cast<rad>((tg->ifield[GETINDEX(i1,i2,i3,n1,n2,n3)])); 
    /*if(rf !=nullptr)
    {
       std::cout<<"res "<<std::endl;
       std::cout<<"rad field temp "<<rf->radtemp<<" "<<std::endl;
    }
    else
      std::cout<<"nullptr"<<std::endl;   
    */


	result=simmesh.m_na;
	return result;
}

int CreateSweepGraph()
{
int result;
  int i1,i2,i3;
  RadInputData data;
  radtmesh simmesh(data);
  simmesh.creategrid();
  simmesh.initgrid();
 

  //std::shared_ptr<grid> tg=simmesh.octants[0]->mgrid[1];
   // std::shared_ptr<rad>rf;
   // std::shared_ptr<rad>tnf;
 int n1=simmesh.nx;//tg->n1;
    int n2=simmesh.ny;//tg->n2;
    int n3=simmesh.nz;//tg->n3; //number of cells in each direction

i1=0;
i2=0;
i3=0;

   // rf= std::static_pointer_cast<rad>((tg->ifield[0]));
    //rf= std::static_pointer_cast<rad>((tg->ifield[GETINDEX(i1,i2,i3,n1,n2,n3)])); 
   /* if(rf !=nullptr)
       std::cout<<rf->radtemp<<" "<<std::endl;
    if(rf!=nullptr && rf->nbr[0]!=nullptr)
    {
      tnf=std::static_pointer_cast<rad>(rf->nbr[0]);
      std::cout<<rf->radtemp<<" "<< tnf->radtemp<<std::endl;
    }
    else
      std::cout<<"nullptr"<<std::endl;*/

      std::cout<<"setup radsweepgraph"<<std::endl;
 simmesh.radsetupsweepgraph();

std::cout<<"test rad pointers"<<std::endl;
  std::shared_ptr<rad>rf;
  std::shared_ptr<rad>tnf;
  std::shared_ptr<rad>tnf2;
  std::shared_ptr<rad>tnf3;
  std::shared_ptr<rad>crad =std::make_shared<rad>(simmesh.m_nf,simmesh.m_na,1);
  std::shared_ptr<rad>crad1 =std::make_shared<rad>(simmesh.m_nf,simmesh.m_na,1);
  std::shared_ptr<rad>crad2 =std::make_shared<rad>(simmesh.m_nf,simmesh.m_na,1);
  std::shared_ptr<rad>crad3 =std::make_shared<rad>(simmesh.m_nf,simmesh.m_na,1);
  std::vector<std::shared_ptr<rad >> rads;
  rf=crad1;
  crad->radtemp=30;
  crad1->radtemp=60;
  crad2->radtemp=70;
  crad3->radtemp=80;

  crad1->nbr[0]=crad1;
  crad1->nbr[1]=crad2;
  crad1->nbr[2]=crad3;
  crad1->nbr[3]=crad1;

  tnf=std::static_pointer_cast<rad>(crad1->nbr[1]);
  tnf2=std::static_pointer_cast<rad>(crad1->nbr[2]);
  tnf3=std::static_pointer_cast<rad>(crad1->nbr[3]);
  //rf=std::static_pointer_cast<rad>(crad->copy());
  rf->radtemp=50;
  rads.push_back(rf);
  rads.push_back(crad);
  std::cout<<"end test rad pointers"<<std::endl;
  std::cout<<"end test rad pointers"<<tnf->radtemp<<std::endl;
  std::cout<<"end test rad pointers"<<tnf2->radtemp<<std::endl;
  std::cout<<"end test rad pointers"<<tnf3->radtemp<<std::endl;



 simmesh.setnbrcellptrs();
 
  result=simmesh.m_na;

  /*std::shared_ptr<radtmesh> rtsimmesh = std::make_shared<radtmesh>(10,1,1);
  radraytsimulation rtsimulation(rtsimmesh);
  rtsimulation.create("");
  rtsimulation.initgrids();
  result=12;*/
 
	return result;
}

int SingleSweepTest()
{
int result;
RadInputData data;

  
	
  data.display();
   std::cout << "test read input file" << std::endl;
   if (data.readInputFile("input.txt")) {
        data.display();  // Display parsed data
    } else {
        std::cerr << "Failed to read input file!" << std::endl;
    }






  radtmesh simmesh(data);
  simmesh.creategrid();
  simmesh.initgrid();
  simmesh.setnbrcellptrs();
        
 simmesh.radsetupsweepgraph();
simmesh.computejnu();
simmesh.computerads();
simmesh.computeradq();
simmesh.computeradtemp();

   simmesh.writeradjfieldvtk( 1000,0);
    simmesh.writeradjfieldvtk( 1000,2);
    simmesh.writeradjfieldvtk( 1000,4);


    simmesh.writeradifieldvtk( 1000,0,0);
    simmesh.writeradifieldvtk( 1000,2,0);
    simmesh.writeradifieldvtk( 1000,4,0);


 std::cout<<"radsweepgraph now setup"<<std::endl;
 simmesh.integratefields();

   simmesh.writeradjfieldvtk( 2000,0);
    simmesh.writeradjfieldvtk( 2000,2);
    simmesh.writeradjfieldvtk( 2000,4);


    simmesh.writeradifieldvtk( 2000,0,0);
    simmesh.writeradifieldvtk( 2000,2,0);
    simmesh.writeradifieldvtk( 2000,4,0);





  result=simmesh.m_na;
	return result;
}



int IntegrateFieldsTest()
{
  int result;
  RadInputData data;
  radtmesh simmesh(data);


  result=simmesh.m_na;
	return result;
}

int AdvanceFieldsTest()
{
  int result;
  RadInputData data;
  radtmesh simmesh(data);


  result=simmesh.m_na;
	return result;
}



int SourceTest()
{
  int result;
  RadInputData data;
  radtmesh simmesh(data);


  result=simmesh.m_na;
	return result;
}

int BoundaryTest()
{
  int result;
  RadInputData data;
  radtmesh simmesh(data);


  result=simmesh.m_na;
	return result;
}


int InstantiateRadtsimulation()
{
	int result;
  RadInputData data;

  
	
  data.display();
   std::cout << "test read input file" << std::endl;
   if (data.readInputFile("input.txt")) {
        data.display();  // Display parsed data
    } else {
        std::cerr << "Failed to read input file!" << std::endl;
    }
  
  std::shared_ptr<radtmesh> rtsimmeshtest;

  //instantiate radiation mesh
  std::shared_ptr<radtmesh> rtsimmesh = std::make_shared<radtmesh>(data);
  radraytsimulation rtsimulation(rtsimmesh);

  std::shared_ptr<mesh> rtmeshbasePtr = std::static_pointer_cast<mesh>(rtsimmesh);

  rtsimulation.create("");
  //radraytsim.initgrids();
  //radraytsim.rungrids();
  //double x =rtsimmesh.getoctpos();
  //cout << "inst simtmesh " << x << '\n';

  rtsimmesh->m_na=24;
  std::cout << "\n\n\nfrom rtsimmesh" << rtsimmesh->m_na << " \n\n\n";

  rtsimmeshtest = std::dynamic_pointer_cast<radtmesh> (rtsimulation.simmesh);
  rtsimmeshtest->m_na=48;
  std::cout << "\n\n\nfrom rtsimmesh" << rtsimmesh->m_na << " \n\n\n";

  //rtsimmeshtest=std::shared_ptr<radtmesh> rtsimulation.simmesh;

  //cout << rtsimmeshtest.m_na << " \n";
	
	//result=rtsimmesh.m_na;
  result=12;
	return result;
}


// Returns n! (the factorial of n).  For negative n, n! is defined to be 1.
int Factorial(int n) {
  int result = 1;
  for (int i = 1; i <= n; i++) {
    result *= i;
  }

  return result;
}

// Returns true if and only if n is a prime number.
bool IsPrime(int n) {
  // Trivial case 1: small numbers
  if (n <= 1) return false;

  // Trivial case 2: even numbers
  if (n % 2 == 0) return n == 2;

  // Now, we have that n is odd and n >= 3.

  // Try to divide n by every odd number i, starting from 3
  for (int i = 3;; i += 2) {
    // We only have to try i up to the square root of n
    if (i > n / i) break;

    // Now, we have i <= n/i < n.
    // If n is divisible by i, n is not prime.
    if (n % i == 0) return false;
  }

  // n has no integer factor in the range (1, n), and thus is prime.
  return true;
}
