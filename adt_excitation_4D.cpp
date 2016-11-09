#include <iostream>
#include <string>
#include <cmath>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <map>
#include <random>

using namespace std;

const double pi = 3.1415926535897; 

void printProgBar( int percent ){
  std::string bar;

  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }

  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%     " << std::flush;
}

float tuneAmp (float amplitude, float a, float b, float qx0){
  float dtune;
  dtune = a + b*amplitude*amplitude;
  return (qx0 + dtune);
}

float oneTurnMapX(float x0, float px0, float y0, float tune, float k_sext, float k_oct){
  float x;
  x = cos(2*pi*tune)*x0 + sin(2*pi*tune)*(px0 + k_sext*x0*x0 - k_oct*x0*x0*x0);
  return (x);
}

float oneTurnMapPx(float x0, float px0, float y0, float tune, float k_sext, float k_oct){
  float px;
  px = -sin(2*pi*tune)*x0 + cos(2*pi*tune)*(px0 + k_sext*(x0*x0 - y0*y0) - k_oct*x0*x0*x0);
  return (px);
}

float oneTurnMapY(float y0, float py0, float x0, float tune, float k_sext, float k_oct){
  float y;
  y = cos(2*pi*tune)*y0 + sin(2*pi*tune)*(py0 - 2.0*k_sext*x0*y0 - k_oct*x0*x0*x0);
  return (y);
}

float oneTurnMapPy(float y0, float py0, float x0, float tune, float k_sext, float k_oct){
  float py;
  py = -sin(2*pi*tune)*y0 + cos(2*pi*tune)*(py0 - 2.0*k_sext*x0*y0 - k_oct*x0*x0*x0);
  return (py);
}

float ADT_kick(float ampl_adt, float tune_adt, float adt_phase, float turn){
  float p;
  p = ampl_adt*cos(2*pi*tune_adt*turn + adt_phase);
  return (p);
}

float ADT_ramp(float tune_adt_0, float tune_adt_f, int Nturns, int turn){
  float tune_adt = tune_adt_0 + (tune_adt_f - tune_adt_0)/Nturns*turn;
  return (tune_adt);
}

float ADT_step_ramp(float tune_adt, float tune_adt_0, float tune_adt_f, int Nturns, int turn, int steps){
  float tune_step = (tune_adt_f - tune_adt_0)/steps;
  int turn_step = Nturns/steps;
  if (turn%turn_step == 0)
    { 
      tune_adt = tune_adt + tune_step;
    }
  return (tune_adt);
}

float ADT_amp_ramp_up(float ampl_adt_f, int ramp_up_turns, int turn){
  float ampl_adt = ampl_adt_f/ramp_up_turns*turn;
  return (ampl_adt);
}

float ADT_amp_ramp_down(float ampl_adt_f, int rampDownStart, int Nturns, int turn){
  float ampl_adt = ampl_adt_f - ampl_adt_f/(Nturns - rampDownStart)*(turn - rampDownStart);
  return (ampl_adt);
}

float ADT_phase(float tune0, float tune, int turn){
  float ADT_phase = (tune-tune0)*turn - tune0;
  return (ADT_phase);
}

int main ()
{
cout << "\n";
  cout << "Active halo control simulations\n";
  cout << "A simple model to prove the principle\n";
  cout << "\n";

  int Npart = 1;
  float time = 0.1; // seconds
  int rev_f = 11.2e3;
  int Nturns = time*rev_f;
  int Nsteps = 1;
  int rampUpEnd = 0.0;
  int rampDown = 0.0;
  float adt_gain = 0.03;
  float adt_max_kick = 0.3; // at 450 GeV
  float ampl_adt_f = 0.684e-2*0.0;
  float ramp_turns = 100.0;
  float tune_adt_0 = 64.290;
  float tune_adt_f = 64.300;
  float qx0 = 64.28;
  float qy0 = 59.31;

  float beam_size = 2.2;
  float collimator_cut = 6.0;

  float k_sext = 0.0;
  float k_oct = -0.93e-2*0.0;
  
  // Detuning with amplitude coefficients measured at 450 GeV
  //float a = 1.7e-3;
  //float b = 0.52e-3;

  cout << "Enter number of turns\n";
  cin >> Nturns;
  cout << "Enter number of particles\n";
  cin >> Npart;
  cout << "Beam size (normalized units)\n";
  cin >> beam_size;

  cout << "Initial ADT frequency (tune units)\n";
  cin >> tune_adt_0;
  cout << "Final ADT frequency (tune units)\n";
  cin >> tune_adt_f;
  cout << "Number of frequency steps\n";
  cin >> Nsteps;

  cout << "ADT ramp up (turns)\n";
  cin >> rampUpEnd;
  cout << "ADT ramp down (turns)\n";
  cin >> rampDown;
  float rampDownStart = Nturns - rampDown;

  cout << "Number of turns = " << Nturns << endl;
  cout << "Number of particles = " << Npart << endl;
  cout << "\n";
  cout << "//// ADT settings //// \n";
  cout << "ADT amplitude = " << ampl_adt_f << endl;
  cout << "Initial ADT tune = " << tune_adt_0 << endl;
  cout << "Final ADT tune = " << tune_adt_f << endl;
  cout << "\n";

  ofstream myfile0, myfile1, myfile2, myfile3, myfile4, myfile5;
  ofstream myfile6, myfile7, myfile8, myfile9, myfile10, myfile11;

  myfile0.open ("dist_at_initial.dat");
  myfile1.open ("dist_at_10.dat");
  myfile2.open ("dist_at_20.dat");
  myfile3.open ("dist_at_30.dat");
  myfile4.open ("dist_at_40.dat");
  myfile5.open ("dist_at_50.dat");
  myfile6.open ("dist_at_60.dat");
  myfile7.open ("dist_at_70.dat");
  myfile8.open ("dist_at_80.dat");
  myfile9.open ("dist_at_90.dat");
  myfile10.open ("dist_at_final.dat");

  ofstream myfile20;
  myfile20.open ("ramp.dat");
  myfile20 << "# [1]Turn " << " [2]ADT_tune " << " [3]ADT_ampl " << " [4]ADT_phase " << endl;
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0,beam_size/sqrt(2.0));

  float x, px, x0, px0, y, py, y0,  py0;
  float amplitudeX, amplitudeY, tune, tune_adt, ampl_adt, adt_phase=0.0, tune_adt_prev;
  bool ramp = true;

  int j = 1;
  while (j <= Npart){    
    int  n = 1;
    double numberX = distribution(generator);
    double numberPX = distribution(generator);    
    double numberY = distribution(generator);
    double numberPY = distribution(generator);
    x0 = numberX;
    px0 = numberPX;
    y0 = numberY;
    py0 = numberPY;
    //p0 = 0.0;
    amplitudeX = sqrt(x0*x0 + px0*px0);
    amplitudeY = sqrt(y0*y0 + py0*py0);
    
    tune_adt = tune_adt_0;
    myfile0 << j << " " << n << " " << x0 << " " << px0 << " " << y0 << " " << py0 << " " << amplitudeX <<" " << amplitudeY << " "<< tune_adt << endl;

    while (n <= Nturns){
      amplitudeX = sqrt(x0*x0 + px0*px0);
      amplitudeY = sqrt(y0*y0 + py0*py0);
      tune = qx0;
      //tune_adt = ADT_ramp(tune_adt_0,tune_adt_f,Nturns,n); // Linear ramp

      tune_adt_prev = tune_adt;
      tune_adt = ADT_step_ramp(tune_adt_prev, tune_adt_0, tune_adt_f, Nturns, n, Nsteps); // Step ramp
      
      if (tune_adt != tune_adt_prev)
	{
	  adt_phase =  ADT_phase(tune_adt_0, tune_adt, n);
	}

      if(n<=rampUpEnd)
	{
	  ampl_adt = ADT_amp_ramp_up(ampl_adt_f, rampUpEnd, n);
	}
      else
	{
	  ampl_adt = ampl_adt_f;
	}
      if(n>=rampDownStart)
	{
	  ampl_adt = ADT_amp_ramp_down(ampl_adt_f, rampDownStart, Nturns, n);
	}

      if (ramp)
	{
	  myfile20 <<  n  << " " << tune_adt  << " " << ampl_adt<<  " " << adt_phase << endl;
	}

      //cout <<  n << " " << (tune_adt_f - tune_adt_0)/Nsteps << endl;

      x = oneTurnMapX(x0,px0,y0,tune, k_sext, k_oct);
      px = oneTurnMapPx(x0,px0,y0,tune, k_sext, k_oct) + ADT_kick(ampl_adt,tune_adt,adt_phase,n);
      y = oneTurnMapY(y0,py0,x0,tune, k_sext, k_oct);
      py = oneTurnMapPy(y0,py0,x0,tune, k_sext, k_oct);


      if (n == 1*Nturns/10){
	myfile1 << j << " " << n << " " << x << " " << px << " " << y << " " << py << " "  << amplitudeX << " " << amplitudeY << " " <<  tune_adt << endl;
      }

      if (n == 2*Nturns/10){
	myfile2 << j << " " << n << " " << x << " " << px << " " << y << " " << py << " "  << amplitudeX << " " << amplitudeY << " " <<  tune_adt << endl;
      }

      if (n == 3*Nturns/10){
	myfile3 << j << " " << n << " " << x << " " << px << " " << y << " " << py << " "  << amplitudeX << " " << amplitudeY << " " <<  tune_adt << endl;
      }

      if (n == 4*Nturns/10){
	myfile4 << j << " " << n << " " << x << " " << px << " " << y << " " << py << " "  << amplitudeX << " " << amplitudeY << " " <<  tune_adt << endl;
      }

      if (n == 5*Nturns/10){
	myfile5 << j << " " << n << " " << x << " " << px << " " << y << " " << py << " "  << amplitudeX << " " << amplitudeY << " " <<  tune_adt << endl;
      }

      if (n == 6*Nturns/10){
	myfile6 << j << " " << n << " " << x << " " << px << " " << y << " " << py << " "  << amplitudeX << " " << amplitudeY << " " <<  tune_adt << endl;
      }

      if (n == 7*Nturns/10){
	myfile7 << j << " " << n << " " << x << " " << px << " " << y << " " << py << " "  << amplitudeX << " " << amplitudeY << " " <<  tune_adt << endl;
      }

      if (n == 8*Nturns/10){
	myfile8 << j << " " << n << " " << x << " " << px << " " << y << " " << py << " "  << amplitudeX << " " << amplitudeY << " " <<  tune_adt << endl;
      }

      if (n == 9*Nturns/10){
	myfile9 << j << " " << n << " " << x << " " << px << " " << y << " " << py << " "  << amplitudeX << " " << amplitudeY << " " <<  tune_adt << endl;
      }
     
      x0 = x;
      px0 = px;
      y0 = y;
      py0 = py;

      n++;
    }  	

    ramp = false;
    myfile20.close();  

    myfile10 << j << " " << n << " " << x0 << " " << px0 << " " << y0 << " " << py0 << " " << amplitudeX <<" " << amplitudeY << " "<< tune_adt << endl;

    printProgBar(1.0*j/Npart*100);

    j++;
  }

  myfile0.close();
  myfile1.close();
  myfile2.close(); 
  myfile3.close();
  myfile4.close();
  myfile5.close(); 
  myfile6.close();
  myfile7.close(); 
  myfile8.close();
  myfile9.close();
  myfile10.close(); 
  return 0;
}

