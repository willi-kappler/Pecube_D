#define HE_NODES 202
#define HE_STEPS 700002

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "aft.h"

using namespace std;





int NSTEPS = 20;
double MILYEARS = 50.0;                 // time padding for thermochronometer work
double TIME_STEP = 20000.0;             // time step size in FEM calculation, units = years. converted to sec in program.
double B = 3.160e-22;                   // material constant (My)
double EBAR = 49.77;                    // acivation energy (kcal / mol)
double R = 0.0019872;                   // gas constant (kcal / (mol * deg))
double DPAR = 1.8;
double ALPHA = -0.12327;
double BETA = -11.988;
double MIN_DETECT = 0.13;
double C[4] = {-19.844, 0.38951, -51.253, -7.6423};






// FUNCTION: Determine The Aft Age Of Each Tracked Node
int get_aft_age_(float* time_orig,float* temp_orig,int* m,float* aftAge)
{
  //  extern int NSTEPS;
  //  extern double MILYEARS;
  // time information must be in reverse order from pecube:
  // for ex. 1.0My, 2.0My, 3.0My, 4.0My, ...
  // for the temperature information this must not be changed

  int i = 0;
  double deltime = 0.0;
  double time_int = 0.0;
  double age = 0.0;
  double* temp_pad;
  double* time_pad;
  double* timein;
  double* tempin;

  timein = new double[*m];
  tempin = new double[*m];
  
  deltime = time_orig[*m - 1];
  for(i = 0;i < *m;i++)
  {
    timein[i] = deltime - time_orig[i];
  }

  // cout << "deltime: " << deltime << endl;
  
  // convert temperature values to Kelvin
  for(int i = 0;i < *m;i++)
  {
    tempin[i] = temp_orig[i] + 273.0;
  }
  
  // pad temp history with the first recored temp of the node deg C values for MILYEARS 
  // in NSTEPS to allow for simulation of unrest ages
  time_pad = new double[*m + NSTEPS];
  temp_pad = new double[*m + NSTEPS];

  time_int = MILYEARS / NSTEPS;

  for(i = 0;i < (NSTEPS - 1);i++)
  {
    time_pad[i] = timein[0] + MILYEARS - (double)(i + 1) * time_int;
  }
  
  for(i = 0;i < *m;i++)
  {
    time_pad[i + NSTEPS - 1] = timein[i];
  }



  for(i = 0;i < (NSTEPS - 1);i++)
  {
    temp_pad[i] = tempin[0];
  }
  
  for(i = 0;i < *m;i++)
  {
    temp_pad[i + NSTEPS - 1] = tempin[i];
  }

  // pass temp history to annealing function
  aft_age_fun(time_pad,temp_pad,(*m + NSTEPS - 1),&age);

  *aftAge = (float) age;
  
  delete [] time_pad;
  delete [] temp_pad;
  delete [] timein;
  
  return (0);
}



// FUNCTION: Determine The Aft Age From A Temperature History With The Methods Of Ketchum Et Al. 2001. Input Temp History With  My bp.
int aft_age_fun(double* time_pad,double* temp_pad,int nin,double *age2)
{
  int nsubstp = nin;
  double age = 0.0;
  double* subtime = new double[nsubstp];
  double* subtemp = new double[nsubstp];

  // divide temperature history into substeps (components) of constant length

  // make substeps arrays
  for(int i = 0;i < nsubstp;i++)
  {
    subtime[i] = time_pad[i];
    subtemp[i] = temp_pad[i];
  }

  anneal(subtemp,subtime,nsubstp,&age);
  *age2 = age;

  delete [] subtime;
  delete [] subtemp;

  return(0);	
}


// FUNCTION: Divide The Temperature History Into Substeps To For Generation Of Components Of Track Annealing
int make_substeps(double* time_pad,double* temp_pad,int nsubstp,double* subtemp,double* subtime)
{
  for(int i = 0;i < nsubstp;i++)
  {
    subtime[i] = time_pad[i];
    subtemp[i] = temp_pad[i];
  }

  return (0);
}



// FUNCTION: Anneal The Mean Track Length For Given Component
int anneal(double* subtemp,double* subtime,int nsubstp,double* age2)
{
  int i = 0;
  int imore = 0;
  int nsubsub = 0;
  double deltime = 0.0;
  double maxdeltime = 0.0;
  double age = 0.0;
  int *addmore;
  //  int *addmore_new;
  int *fst_subsub;
  double *sstime;
  double *sstemp;
  double *sstime_new;
  double *sstemp_new;
  double *t_length;

  // make subsub-steps and ensure that the delta temp is within the criteria.
  // determine the number of subsub-steps between each substep

  fst_subsub = new int[nsubstp];

  // allocate space for and make first subsub-step intervals based on max change in time
  maxdeltime = 0.001;

  // make sure maxdeltime is less than the time interval of substeps
  for(i = 0;i < (nsubstp - 1);i++)
  {
    if (subtime[i] - subtime[i + 1] < 0.0) {
      fprintf(stderr, "\n in function anneal: invalid time configuration in input: \n");
      fprintf(stderr, "( %16.8f - %16.8f ) < 0.0, at index: %d\n", subtime[i], subtime[i + 1], i);
      *age2 = 0.0;
      return 0;
    }

    while(maxdeltime >= (subtime[i] - subtime[i + 1]))
    {
      maxdeltime = 0.7 * maxdeltime;
      fprintf(stderr,"\nmaxdeltime is > substep length\n newmaxdeltime=%f",maxdeltime);

      // fprintf(stderr,"\nmaxdeltime is > substep length\n newmaxdeltime=%f",maxdeltime,subtime[i] - subtime[i + 1]);
    }	
  }		

  deltime = subtime[0] - subtime[nsubstp - 1];
  nsubsub = static_cast<int>(floor(deltime / maxdeltime) + nsubstp);

  sstime = new double[nsubsub];
  sstemp = new double[nsubsub];
  addmore = new int[nsubsub];

  sstemp[0] = subtemp[0];
  sstime[0] = subtime[0];

  // check that deltemp between each subsub-step is less than 8c everywhere and 3.5c when
  // within 10c of the total Annealing temperature using the subsub-step to define 
  // a cooling rate. if conditions not met, add more substeps.
  imore = 1;

  while(imore > 0)
  {
    check_tot_anneal(&imore,addmore,sstemp,sstime,nsubsub,subtemp,subtime,nsubstp,fst_subsub);

    // if additional subsub-steps are needed (imore>0), reallocate space for time and temp arrays 
    // and cut subsub-step in half at the problem step
    if(imore > 0)
    {
      nsubsub = nsubsub + imore;
      
      sstime_new = new double[nsubsub];
      sstemp_new = new double[nsubsub];
      // addmore_new = new int[nsubsub];		
      
      add_subsub_steps(nsubsub,imore,sstime_new,sstemp_new,sstemp,sstime,addmore);
      
      delete [] sstime;
      delete [] sstemp;
      delete [] addmore;
      
      sstime = new double[nsubsub];
      sstemp = new double[nsubsub];
      addmore = new int [nsubsub];
      
      for(i = 0;i < nsubsub;i++)
      {
        sstime[i] = sstime_new[i];
        sstemp[i] = sstemp_new[i];
      }	
      
      delete [] sstime_new;
      delete [] sstemp_new;
    }	
  }	

  // component substeps are subtime,subtemp with nsubstp components. Subsub-steps for
  // numerical integration are sstime and sstemp where the ith substep begins with subsub-step
  // fst_subsub[0][i].  check that thses beginings are correct

  // allocate normalized track length array
  t_length = new double[nsubstp];

  // anneal all components through time
  anneal_comp(t_length,nsubstp,fst_subsub,nsubsub,sstime,sstemp);

  // calculate the age of the sample
  calc_age_(t_length,&nsubstp,subtime,&age);
  *age2 = age;	

  delete [] fst_subsub;
  delete [] addmore;
  delete [] t_length;
  delete [] sstime;
  delete [] sstemp;

  return (0);
}



// FUNCTION: Calculate The Age Of The Sample Using Track Density
int calc_age_(double* t_length,int* nsubstp,double* subtime,double* age)
{
  int icomp = 0;
  double rhocomp = 0.0;
  double deltime = 0.0;
  // double wt_decay = 0.0;
  double rhost = 0.893;
  // double u_decay = 0.000155125;

  *age = 0.0;

  for(icomp = 0;icomp < *nsubstp - 1;icomp++)
  {
    if(t_length[icomp] >= 0.765)
    {
      rhocomp = 1.600 * t_length[icomp] - .600;
    }
    else if((t_length[icomp] <= 0.0) || (t_length[icomp] != t_length[icomp]))
    {
      rhocomp = 0.0;
    }	
    else
    {
      rhocomp = 9.205 * pow(t_length[icomp],2.0) - 9.157 * t_length[icomp] + 2.269;
    }
    deltime = fabs(subtime[icomp + 1] - subtime[icomp]);
    // wt_decay = (exp(subtime[icomp] * u_decay) - exp(subtime[icomp + 1] * u_decay)) / u_decay;

    *age = *age + rhocomp * deltime;
  }

  //  cout << endl << "calc_age_: age: " << (*age) << ", rhocomp: " << rhocomp << ", deltime: " << deltime << endl;

  *age = *age / (rhost);

  return (0);
}



// FUNCTION:  Anneal All Components Through Time
int anneal_comp(double* t_length,int nsubstp,int* fst_subsub,int nsubsub,double* sstime,double* sstemp)
{
  int isubsub = 0;
  int icomp = 0;
  int ipoint = 0;
  int icatch = 0;
  int i = 0;
  double rhs = 0.0;
  double deltime = 0.0;
  double prev_length = 0.0;
  double rmro = 0.0;
  double kappa = 0.0;
  double temp = 0.0;
  double lhs = 0.0;
  double time_eq = 0.0;
  double time = 0.0;
  //  double test = 0.0;
  double mpoint = 0.0;

  // define constants for annealing equation
  // extern double DPAR;
  // extern double ALPHA;
  // extern double BETA;
  // extern double MIN_DETECT;
  // extern double C[4];

  fst_subsub[nsubstp - 1] = nsubsub - 1;	// add this value for the loop below to find the start point
                                          // of integration for the last component

  for(icomp = 0;icomp < nsubstp - 1;icomp++)
  {
    t_length[icomp] = 1.0;
    prev_length = 1.0;
    
    // pick the midpoint of the comp timestep as the beginning of integration
    mpoint = sstime[fst_subsub[icomp]] - .5 * (fabs(sstime[fst_subsub[icomp]] - sstime[fst_subsub[icomp + 1]]));
    icatch = 1;
    for(i = fst_subsub[icomp];i < fst_subsub[icomp + 1];i++)
    {
      if(sstime[i] < mpoint)
      {
        ipoint = i - 1;
        icatch = 0;
        break;
      }
    }	
    if(icatch == 1)
    {
      fprintf(stderr,"\n\n#################  ERROR: Did Not Find Midpoint At Which To Begin Integration  #################\n\n");
      ipoint = fst_subsub[icomp];
    }	
    for(isubsub = ipoint;isubsub < nsubsub - 1;isubsub++)
    {
      // define temp as the average over the subsubstep
      temp = sstemp[isubsub] - .5 * (sstemp[isubsub] - sstemp[isubsub + 1]);
      
      // Find the equivalent time to reach the length at the begining of the time step at temperature = temp
      lhs = (pow(((1.0 - pow(prev_length,BETA)) / BETA),ALPHA) - 1.0) / ALPHA;
      time_eq = exp(((lhs - C[0]) / C[1]) * (log(1.0 / temp) - C[3]) + C[2]);
      deltime = fabs(sstime[isubsub] - sstime[isubsub + 1]) * 3.15578e13;
      time = deltime + time_eq;
      
      // solve for the new length with time
      rhs = C[0] + C[1] * ((log(time) - C[2]) / (log(1.0 / temp) - C[3]));
      t_length[icomp] = pow(1.0 - BETA * pow(ALPHA * rhs + 1.0,1 / ALPHA),1 / BETA);
      prev_length = t_length[icomp];

      // Having problem with the track lengths getting so small that the rhs(?)
      // values are bigger than the machine can handle. So, if the tracks get really
      // small, set the length to zero and break out of the loop. There is some
      // concern that the length cutoff of .01 might not be small enough when
      // converted for differnt compositions, but I doubt it.
      if((t_length[icomp] < .01) || (t_length[icomp] != t_length[icomp]))
      {
        t_length[icomp] = 0.0;
        break;
      }	
    }	
    // convert track lengths to track length for an apatite of composition dpar
    rmro = 1.0 - exp(.647 * (DPAR - 1.75) - 1.834);
    kappa = 1.0 - rmro;
    
    // check that length is > min length, also that converted length is > 0, ie can be taken to the power kappa
    if(pow((t_length[icomp] - rmro) / (1.0 - rmro),kappa) < MIN_DETECT)
    {
      t_length[icomp] = 0.0;
    }
    if((t_length[icomp] - rmro) / (1.0 - rmro) < 0.0)
    {
      t_length[icomp] = 0.0;
    }	
    else
    {
      t_length[icomp] = pow((t_length[icomp] - rmro)/(1.0 - rmro),kappa);
    }	

    // for isothermal run can check the integrated soln against the analytical soln
    rhs = C[0] + C[1] * ((log(3.15578e13 * fabs(sstime[fst_subsub[icomp]] - sstime[nsubsub - 1])) - C[2]) / (log(1.0 / 273.0) - C[3]));
    // test = pow(1.0 - BETA * pow(ALPHA * rhs + 1.0,1 / ALPHA),1 / BETA);

  }

  return (0);
}



// FUNCTION: Add Additional Subsub Steps Where The Del Temp Criteria Was Not Met. Add One Extra Time Step.
int add_subsub_steps(int nsubsub,int imore,double* sstime_new,double* sstemp_new,double* sstemp,double* sstime,int* addmore)
{
  int i = 0;
  int isubsub = 0;
  double dtemp = 0.0;
  double dtime = 0.0;

  for(i = 0;i < (nsubsub - imore);i++)
  {
    isubsub = 0;
    sstime_new[i + isubsub] = sstime[i];
    sstemp_new[i + isubsub] = sstemp[i];
    if(i == addmore[isubsub])
    {
      dtemp = sstemp_new[i + isubsub - 1] - sstemp_new[i + isubsub];
      dtime = sstime_new[i + isubsub - 1] - sstime_new[i + isubsub];
      sstemp_new[i + isubsub + 1] = sstemp_new[i + isubsub];
      sstime_new[i + isubsub + 1] = sstime_new[i + isubsub];
      sstemp_new[i + isubsub] = sstemp_new[i + isubsub - 1] + .5 * dtemp;
      sstime_new[i + isubsub] = sstime_new[i + isubsub - 1] + .5 * dtime;
      isubsub = isubsub + 1;
    }	
  }

  return (0);
}



// FUNCTION: Check That The Deltemp Between Each Subsub Step Is Less Than 
// 8C Everywhere And 3.5C When Within 10C Of The Total Annealing Temperature
// Using The Subsub Step To Define A Cooling Rate. 
int check_tot_anneal(int* imore,int* addmore,double* sstemp, double* sstime,int nsubsub,double* subtemp,double* subtime,int nsubstp,int* fst_subsub)
{
  int i;
  int isubstp = 0;
  int index = 0;
  int index2 = 0;
  int insert = 0;
  double tot_anneal = 0.0;
  double dtemp = 0.0;
  double new_time = 0.0;
  double deltemp = 0.0;
  double deltime = 0.0;

  *imore = 0;

  // make initial array of subsub-step temps and times inserting extra subsub-steps on both
  // sides of the substep boundaries to ensure that a subsub-step will have the same values as
  // the substep and that the maximum delta time is not exceeded
  for(i = 0;i < nsubsub;i++)
  {
    addmore[i] = 0;
  }	

  for(i = 0;i < nsubsub - nsubstp;i++)
  {
    new_time = subtime[0] - (double)i * (subtime[0] - subtime[nsubstp - 1]) / (double)(nsubsub - nsubstp);
    index2 = 0;
    insert = 0;
    while(new_time <= subtime[isubstp + 1])
    {
      index2 = index2 + 1;
      isubstp = isubstp + 1;
      insert = 1;
      fst_subsub[isubstp] = index;  // set starting subsubstep for isubstp substp
      if(index2 > 1) fprintf(stderr,"\n\n################# ERROR: Shouldn't Skip More Than One SubStep #################\n\n");
    }	
    
    deltemp = subtemp[isubstp] - subtemp[isubstp + 1];
    deltime = subtime[isubstp] - subtime[isubstp + 1];
    if(insert == 0)
    {
      sstemp[index] = subtemp[isubstp] - (subtime[isubstp] - new_time) * deltemp / deltime;
      sstime[index] = new_time;
    }
    else if(new_time == subtime[isubstp])
    {
      sstemp[index + 1] = subtemp[isubstp];
      sstime[index + 1] = new_time;
      sstemp[index] = sstemp[index - 1] + .5 * (sstemp[index + 1] - sstemp[index - 1]);
      sstime[index] = sstime[index - 1] + .5 * (sstime[index + 1] - sstime[index - 1]);
      index = index + 1;
      insert = 0;
    }	
    else if(insert == 1) {
      sstime[index + 1] = new_time;
      sstime[index] = subtime[isubstp];
      sstemp[index + 1] = subtemp[isubstp] - (subtime[isubstp] - new_time) * deltemp / deltime;
      sstemp[index] = subtemp[isubstp];
      index = index + 1;
    }	
    
    index = index + 1;
  }

  fst_subsub[0] = 0;

  // add in 2nd to last and last subsub-step
  sstime[nsubsub - 1] = subtime[nsubstp - 1];
  sstemp[nsubsub - 1] = subtemp[nsubstp - 1];
  sstime[nsubsub - 2] = sstime[nsubsub - 3] - .5 * (sstime[nsubsub - 3] - sstime[nsubsub - 1]);
  sstemp[nsubsub - 2] = sstemp[nsubsub - 3] - .5 * (sstemp[nsubsub - 3] - sstemp[nsubsub - 1]);

  if(index != nsubsub - 2)
  {
    fprintf(stderr,"\n\n#################  ERROR: In Making SubSubSteps  #################\n\n");
  }	

  // check that all delta temps are within maximum temp change criteria
  for(i = 1; i < nsubsub; i++)
  {
    tot_anneal = 377.67 * pow(fabs(sstemp[i - 1] - sstemp[i]) / fabs(sstime[i - 1] - sstime[i]),0.019837);
    dtemp = sstemp[i - 1] - sstemp[i];
    if(fabs(dtemp) >= 8.0)
    {
      *imore = *imore + 1;
      addmore[*imore] = i;
    }	
    else if(sstemp[i] <= tot_anneal + 10.0)
    {
      if(sstemp[i] >= tot_anneal - 10)
      {
        if(dtemp >= 3.5)
        {
          *imore = *imore + 1;
          addmore[*imore] = i;
        }
      }
    }	
  }
    
  return (0);
}
