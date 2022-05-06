#ifndef AFT_H
#define AFT_H


extern "C" {

// function prototypes
int get_aft_age_(float* timein,float* tempin,int* m,float* aftAge);
int aft_age_fun(double* time_pad,double* temp_pad,int nin, double* age2);
int make_substeps(double* time_pad,double* temp_pad,int nsubstp,double* subtemp,double* subtime);
int anneal(double* subtemp,double* subtime,int nsubstp,double* age2);
int calc_age_(double* t_length,int* nsubstp,double* subtime,double* age);
int anneal_comp(double* t_length,int nsubstp,int* fst_subsub,int nsubsub,double* sstime,double* sstemp);
int add_subsub_steps(int nsubsub,int imore,double* sstime_new,double* sstemp_new,double* sstemp,double* sstime,int* addmore);
int check_tot_anneal(int* imore,int* addmore,double* sstemp, double* sstime,int nsubsub,double* subtemp,double* subtime,int nsubstp,int* fst_subsub);

}


#endif
