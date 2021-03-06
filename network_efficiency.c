/*
 *  Copyright (C) 2012 Nickolas Fotopoulos, Leo Singer, Alexander Dietz
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <glob.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <lal/DetResponse.h>
#include <lal/LALSimNoise.h>
#include <lal/LALConstants.h>
#include <lal/LALDetectors.h>
#include <time.h>

int debug=1;

ssize_t str2network(LALDetector network[], double reach[], char *str);
static double horizon_distance(double m1, double m2, double f_low, double f_high, double snr_thresh, char detector);
static double rchisq_2(gsl_rng *rng, double lambda);
void Ssq(double *S2, gsl_rng *rng, double beam_fac, double *response, ssize_t network_size);
static void wave_Consts_Test(double *rho2, double *X2, gsl_rng *rng, double lambda);

void compute_efficiency(double *efficiency, LALDetector *network, size_t network_size, double* BNS_horizon_dist, double beam_fac, double m1_min, double m1_max, double m2_min, double m2_max, double rho_thresh, double *distances, size_t dist_bins, size_t samples){
  const static double twopi=2 * M_PI;

  #pragma omp parallel
  {
    /*Create and seed random number generator*/
    gsl_rng *rng=gsl_rng_alloc(gsl_rng_mt19937);
    #ifdef _OPENMP
    gsl_rng_set(rng, omp_get_thread_num());//Seed with thread number.
    #else
    gsl_rng_set(rng, 0);
    #endif

    double *response=malloc(2 * network_size * sizeof(double));

    /*Total detections for each distance step for the given thread.*/
    long *threadTotals=calloc(dist_bins, sizeof(long));

    /*Preallocate variables needed inside loops*/
    size_t j, k, l;
    /*size_t objects are inherently the size of a pointer on the system.*/
    
    #pragma omp for schedule(static)
    for(k=1; k<=samples; k++){
      /* draw orientations and determine network response */
      const double lon=twopi*gsl_rng_uniform(rng);
      const double lat=M_PI_2 - acos(2*gsl_rng_uniform(rng) - 1);
      const double zeta=twopi*gsl_rng_uniform(rng);
      const double BNSchirp=pow( pow(1.4*1.4,3.0)/(1.4+1.4), 0.2);
      const double m1=m1_min+gsl_rng_uniform(rng)*(m1_max-m1_min);
      const double m2=m2_min+gsl_rng_uniform(rng)*(m2_max-m2_min);
      const double mass_correction=pow( pow(m1*m2,3.0)/(m1+m2), 0.2)/BNSchirp;

      for(l=network_size; l--;)
	XLALComputeDetAMResponse(response+2*l, response+2*l+1, network[l].response, lon, lat, zeta, 0.0);

      double S2[network_size];
      Ssq(S2, rng, beam_fac, response, network_size);
      double lambda_const[network_size];
      for(l=network_size; l--;)
	lambda_const[l]=(rho_thresh*rho_thresh-2) * BNS_horizon_dist[l]*BNS_horizon_dist[l] * S2[l] * mass_correction;

      for(j=dist_bins; j--;){
	int successes=0;
	for(l=network_size; l--;){
	  /* r^2=S^2*(D_hor*sqrt(rho_t^2-2)/D)^2 */
          const double lambda=(1/distances[j]) * (1/distances[j]) * lambda_const[l];
          successes+=rchisq_2(rng, lambda)>rho_thresh*rho_thresh;
	}
	threadTotals[j]+=(successes>=2);
      }
    }

    #pragma omp critical
    for(j=dist_bins; j--;)
      efficiency[j]+=threadTotals[j]/((double)samples);

    gsl_rng_free(rng);
    free(threadTotals);
    free(response);
  }
}

int main(int argc, char *argv[]){
  const clock_t beginTimeStamp=clock();
  const static double twopi=2 * M_PI;
  LALDetector network[LAL_NUM_DETECTORS];
  double Dhor[LAL_NUM_DETECTORS];
  double thresh_snr;
  ssize_t network_size;
  double beam_fac=-1.;
  size_t N;
  size_t D_max;
  size_t D_steps;

  /*Parse cmdline arguments*/
  if (argc != 7) {
    fprintf(stderr, "Usage: snr_thresh network Ntrials jet_semiangle_deg max_distance_Mpc distance_steps\n");
    exit(2);
  }

  /*Initilize network and other cmdline arguments*/
  memset(network, 0, LAL_NUM_DETECTORS);
  thresh_snr=strtod(argv[1], NULL);
  network_size=str2network(network, Dhor, argv[2]);
  if (network_size == -1) exit(2);
  N=strtol(argv[3], NULL, 10);
  beam_fac=1 - cos(strtod(argv[4], NULL) * M_PI / 180.);

  /*Set up/initialize the distance parameters*/
  D_max=strtol(argv[5], NULL, 10);
  D_steps=strtol(argv[6], NULL, 10);
  const double d_step=((double)(D_max-1))/(D_steps-1);
  double distances[D_steps];
  for(size_t i=0; i<D_steps; i++)
    distances[i]=1+d_step*i;

  /*Initialize final storage array of efficiency for each distance step*/
  double *efficiency=calloc(D_steps, sizeof(double));

  compute_efficiency(efficiency, network, network_size, Dhor, beam_fac, 1.31, 1.39, 1.31, 1.39, 8, distances, D_steps, N);
  
  

  /*
  //START THREAD HERE
  #pragma omp parallel
  {
    //Create and seed random number generator
    gsl_rng *rng=gsl_rng_alloc(gsl_rng_mt19937);
    #ifdef _OPENMP
    gsl_rng_set(rng, omp_get_thread_num());//Seed with thread number.
    #else
    gsl_rng_set(rng, 0);
    #endif

    double *response=malloc(2 * network_size * sizeof(double));

    //Total detections for each distance step for the given thread.
    long *threadTotals=calloc(D_steps, sizeof(long));

    //Preallocate variables needed inside loops
    size_t j, k, l;
    //size_t objects are inherently the size of a pointer on the system.
    
    #pragma omp for schedule(static)
    for(k=1; k<=N; k++){
      //draw orientations and determine network response
      const double lon=twopi*gsl_rng_uniform(rng);
      const double lat=M_PI_2 - acos(2*gsl_rng_uniform(rng) - 1);
      const double zeta=twopi*gsl_rng_uniform(rng);

      for(l=network_size; l--;)
	XLALComputeDetAMResponse(response+2*l, response+2*l+1, network[l].response, lon, lat, zeta, 0.0);

      double S2[network_size];
      Ssq(S2, rng, beam_fac, response, network_size);

      for(j=D_steps; j--;){
	const double d=1.0+j*d_step;
	int successes=0;
	for (l=network_size; l--;) {
	  //r^2=S^2*(D_hor*sqrt(rho_t^2-2)/D)^2
          const double lambda=(thresh_snr*thresh_snr-2)*(Dhor[l]/d)*(Dhor[l]/d)*S2[l];
          successes+=rchisq_2(rng, lambda)>thresh_snr*thresh_snr;
	}
	threadTotals[j]+=(successes>=2);
      }
    }

    #pragma omp critical
    for(j=D_steps; j--;)
      efficiency[j]+=threadTotals[j]/((double)N);

    gsl_rng_free(rng);
    free(threadTotals);
    free(response);
  }*/

  if(debug){
    size_t i;
    printf("DISTANCE  EFFICIENCY\n");
    for(i=0; i<D_steps; i++)
      printf("%09.5f %09.5f\n", 1.0+i*d_step, 100*efficiency[i]);
  }
  
  free(efficiency);

  return 0;
}

/////####//#####//#///#/////
/////#///////#////##//#/////
/////###/////#////#/#/#/////
/////#///////#////#//##/////
/////#/////#####//#///#/////

static double rchisq_2(gsl_rng *rng, double lambda) {
  /*Generates a random value from a degree 2 
    chi_squared with non-centrality parameter of lambda*/
  double a=sqrt(0.5*lambda);
  const double temp=gsl_ran_gaussian(rng, 1.0)+a;
  const double temp2=gsl_ran_gaussian(rng, 1.0)+a;
  return (temp*temp)+(temp2*temp2);
}

void Ssq(double *S2, gsl_rng *rng, double beam_fac, double *response, ssize_t network_size){
  /* Calculates the antenna factor for each detector for a random source orientation*/
  const double cosiota=1-beam_fac*gsl_rng_uniform(rng);//beam_fac determines max iota.
  const double cosiotasq=cosiota*cosiota;
  const double iotafac=0.25*(1+cosiotasq)*(1+cosiotasq);

  size_t l=network_size;
  for(;l--;) {
    double fplus=response[2*l], fcross=response[2*l+1];
    S2[l]=(fplus*fplus)*iotafac+(fcross*fcross)*cosiotasq;
  }
}

ssize_t str2network(LALDetector network[LAL_NUM_DETECTORS], double reach[LAL_NUM_DETECTORS], char *str) {
  /*Convert string like "HLVK" to an array of LALDetectors.
  Return the size of the network.*/
  size_t k=0;
  while (k < LAL_NUM_DETECTORS && str[k]) {
    /* WARNING masses, frequency cutoffs, and SNR threshold hardcoded! */
    reach[k]=horizon_distance(1.4, 1.4, 40., 1570., 8., str[k]);
    /* fprintf(stderr, "Detector '%c' horizon distance: %g\n", str[k], reach[k]); */
    if (str[k]=='H') {
      network[k++]=lalCachedDetectors[LAL_LHO_4K_DETECTOR];
    }
    else if (str[k]=='L') {
      network[k++]=lalCachedDetectors[LAL_LLO_4K_DETECTOR];
    }
    else if (str[k]=='V') {
      network[k++]=lalCachedDetectors[LAL_VIRGO_DETECTOR];
    }
    else if (str[k]=='K') {
      /* numbers from private communication with Koji Arai */
      LALDetector *detector=network+k;
      LALFrDetector *frDetector=&(detector->frDetector);
      strncpy(frDetector->name, "KAGRA", LALNameLength);
      strncpy(frDetector->prefix, "K1", 3);
      frDetector->vertexLatitudeRadians=2.396511595913414;
      frDetector->vertexLongitudeRadians=0.6354743806511354;
      frDetector->vertexElevation=372.0;
      frDetector->xArmAltitudeRadians=0.0;
      frDetector->xArmAzimuthRadians=1.076693615555302;
      frDetector->yArmAltitudeRadians=0.0;
      frDetector->yArmAzimuthRadians=5.789082595939991;
      detector=XLALCreateDetector(network+k, frDetector, LALDETECTORTYPE_IFODIFF);
      if (!detector) {
        fprintf(stderr, "Failed to create KAGRA detector\n");
        return -1;
      }
      k++;
    }
    else if (str[k] == 'I') {
      /* numbers from Schutz 2011 network FOMs */
      LALDetector *detector=network+k;
      LALFrDetector *frDetector=&(detector->frDetector);
      strncpy(frDetector->name, "Indigo", LALNameLength);
      strncpy(frDetector->prefix, "I1", 3);
      frDetector->vertexLatitudeRadians=1.3098647554849334;
      frDetector->vertexLongitudeRadians=0.33329486135237268;
      frDetector->vertexElevation=0.0;
      frDetector->xArmAltitudeRadians=0.0;
      frDetector->xArmAzimuthRadians=3.9269908169872414;
      frDetector->yArmAltitudeRadians=0.0;
      frDetector->yArmAzimuthRadians=5.497787143782138;
      detector=XLALCreateDetector(network+k, frDetector, LALDETECTORTYPE_IFODIFF);
      if (!detector) {
        fprintf(stderr, "Failed to create Indigo detector\n");
        return -1;
      }
      k++;
    }
    else {
      fprintf(stderr, "Unrecognized site: %c\n", str[k]);
      return -1;
    }
  }
  return k;
}


static double horizon_distance(double m1, double m2, double f_low, double f_high, double snr_thresh, char detector){
  /* Compute horizon distance for a given detector.  After equation (2) from LLOID paper. */
  double (*S)(double);
  double Dhor, f;

  /* Chirp mass in seconds */
  const double Mchirp=LAL_MTSUN_SI*pow(pow(m1*m2, 3)/(m1+m2), 0.2);

  /* Integration step in Hz */
  static const double df=0.02;

  /* Select the function to use for sampling the PSD */
  switch (detector){
    case 'H':
    case 'L':
    case 'I':
      S=XLALSimNoisePSDaLIGOZeroDetHighPower;
    break;
    case 'V':
      S=XLALSimNoisePSDAdvVirgo;
      break;
    case 'K':
      S=XLALSimNoisePSDKAGRA;
      break;
    default:
      fprintf(stderr, "This line should never be reached.\n");
      exit(EXIT_FAILURE);
  }

  for (Dhor=0, f=f_low; f<f_high; f+=df)
    Dhor+=pow(f, (-7.0/3))/S(f)*df;

  Dhor=sqrt((5./6)*Dhor);
  Dhor*=pow(Mchirp, (5./6))*pow(M_PI, (-2./3))*LAL_C_SI/snr_thresh;

  /* Convert from m to Mpc */
  Dhor/=1e6*LAL_PC_SI;

  return Dhor;
}

static void wave_Consts_Test(double *rho2, double *X2, gsl_rng *rng, double lambda){
  /*Takes two pointers, rho2 and X2, and sets the values to random
    rho2: X^2, df=2, ncp=lambda
    X2: X^2, df=2p-2, ncp=0 
    For use as wave consistancy test from Creighton/Anderson
    Sec 7.8.1.4  */

  const static int p=16;
  const double sigma=1;
  const double mu=sigma*sqrt(lambda/(2*p*p));
  double a[p]; double b[p];
  double asum=0; double bsum=0;
  *X2=0;
  double err=2*0.03;

  for (int i=0; i<p; i++){
    a[i]=gsl_ran_gaussian(rng, sigma*(1+err*mu)/sqrt(p))+mu;
    b[i]=gsl_ran_gaussian(rng, sigma*(1+err*mu)/sqrt(p))+mu;
    asum+=a[i]; bsum+=b[i];
  }

  for (int i=0; i<p; i++)
    *X2+=(pow(asum/p-a[i],2) + pow(bsum/p-b[i],2)) * p/pow(sigma,2);

  *rho2=asum*asum+bsum*bsum;
}
