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

int debug=1;

ssize_t str2network(LALDetector[],double[],char*);
static double horizon_distance(double, double, double, double, double, char);
static double rchisq_2(gsl_rng*, double);
void Ssq(double*, gsl_rng*, double, double*, ssize_t);

int main(int argc, char *argv[]) {
  const static double twopi = 2 * M_PI;
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
  thresh_snr = strtod(argv[1], NULL);
  network_size = str2network(network, Dhor, argv[2]);
  if (network_size == -1) exit(2);
  N = strtol(argv[3], NULL, 10);
  beam_fac = 1 - cos(strtod(argv[4], NULL) * M_PI / 180.);

  /*Set up/initialize the distance parameters*/
  D_max=strtol(argv[5], NULL, 10);
  D_steps=strtol(argv[6], NULL, 10);
  double d_step=(D_max-1)/(D_steps-1);

  /*Initialize final storage array of efficiency for each distance step*/
  double *efficiency = calloc(D_steps, sizeof(double));

  //START THREAD HERE
  /*Create and seed random number generator*/
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng, 0);
  double *response = malloc(2 * network_size * sizeof(double));

  /*Total detections for each distance step for the given thread.*/
  long *threadTotals = calloc(D_steps, sizeof(long));

  /*Preallocate variables needed inside loops*/
  size_t j, k, l;
  double S2[network_size];
  /*size_t objects are inherently the size of a pointer on the system.*/
  for (k=1; k<=N; k++) {
    /* draw orientations and determine network response */
    const double lon=twopi*gsl_rng_uniform(rng);
    const double lat=M_PI_2-acos(2*gsl_rng_uniform(rng)-1);
    const double zeta=twopi*gsl_rng_uniform(rng);

    for (l = network_size; l--;)
      XLALComputeDetAMResponse(response + 2 * l, response + 2 * l + 1, network[l].response, lon, lat, zeta, 0.);

    Ssq(S2, rng, beam_fac, response, network_size);

    for(j=D_steps;j--;){
      const double d = 1. + j * d_step;
      int successes = 0;
      for(l=network_size;l--;) {
          const double lambda = (64-2)*(Dhor[l]/d)*(Dhor[l]/d)*S2[l];
          successes += rchisq_2(rng, lambda) > 64;
      }
      threadTotals[j]+=successes>=2;
    }
  }

  printf("This is before the seg fault.\n");
  gsl_rng_free(rng);
  printf("This is after the seg fault.\n");
  free(threadTotals);
  free(response);
  //END THREAD HERE

  free(efficiency);

  return 0;
}

/////####//#####//#///#/////
/////#///////#////##//#/////
/////###/////#////#/#/#/////
/////#///////#////#//##/////
/////#/////#####//#///#/////

static double rchisq_2(gsl_rng *rng, double lambda) {
  double a = sqrt(0.5 * lambda);
  const double temp = gsl_ran_gaussian(rng, 1.0) + a;
  const double temp2 = gsl_ran_gaussian(rng, 1.0) + a;
  return temp * temp + temp2 * temp2;
}

void Ssq(double *S2, gsl_rng *rng, double beam_fac, double *response, ssize_t network_size){
  /* Calculates the antenna factor for each detector for a random source orientation*/
  const double cosiota=1-beam_fac*gsl_rng_uniform(rng);//beam_fac determines max iota.
  const double cosiotasq=cosiota*cosiota;
  const double iotafac=0.25*(1+cosiotasq)*(1+cosiotasq);

  size_t l=network_size;
  for(;l--;) {
    double fplus=response[2*l], fcross=response[2*l+1];
    S2[l]=fplus * fplus * iotafac + fcross * fcross * cosiotasq; //S^2
  }
}

/*
 * Convert string like "HLVK" to an array of LALDetectors.
 * Return the size of the network.
 */
ssize_t str2network(LALDetector network[LAL_NUM_DETECTORS], double reach[LAL_NUM_DETECTORS], char *str) {
  size_t k = 0;
  while (k < LAL_NUM_DETECTORS && str[k]) {
    /* WARNING masses, frequency cutoffs, and SNR threshold hardcoded! */
    reach[k] = horizon_distance(1.4, 1.4, 40., 1570., 8., str[k]);
    /* fprintf(stderr, "Detector '%c' horizon distance: %g\n", str[k], reach[k]); */
    if (str[k] == 'H') {
      network[k++] = lalCachedDetectors[LAL_LHO_4K_DETECTOR];
    }
    else if (str[k] == 'L') {
      network[k++] = lalCachedDetectors[LAL_LLO_4K_DETECTOR];
    }
    else if (str[k] == 'V') {
      network[k++] = lalCachedDetectors[LAL_VIRGO_DETECTOR];
    }
    else if (str[k] == 'K') {
      /* numbers from private communication with Koji Arai */
      LALDetector *detector = network + k;
      LALFrDetector *frDetector = &(detector->frDetector);
      strncpy(frDetector->name, "KAGRA", LALNameLength);
      strncpy(frDetector->prefix, "K1", 3);
      frDetector->vertexLatitudeRadians = 2.396511595913414;
      frDetector->vertexLongitudeRadians = 0.6354743806511354;
      frDetector->vertexElevation = 372.0;
      frDetector->xArmAltitudeRadians = 0.;
      frDetector->xArmAzimuthRadians = 1.076693615555302;
      frDetector->yArmAltitudeRadians = 0.;
      frDetector->yArmAzimuthRadians = 5.789082595939991;
      detector = XLALCreateDetector(network + k, frDetector, LALDETECTORTYPE_IFODIFF);
      if (!detector) {
        fprintf(stderr, "Failed to create KAGRA detector\n");
        return -1;
      }
      k++;
    }
    else if (str[k] == 'I') {
      /* numbers from Schutz 2011 network FOMs */
      LALDetector *detector = network + k;
      LALFrDetector *frDetector = &(detector->frDetector);
      strncpy(frDetector->name, "Indigo", LALNameLength);
      strncpy(frDetector->prefix, "I1", 3);
      frDetector->vertexLatitudeRadians = 1.3098647554849334;
      frDetector->vertexLongitudeRadians = 0.33329486135237268;
      frDetector->vertexElevation = 0.0;
      frDetector->xArmAltitudeRadians = 0.;
      frDetector->xArmAzimuthRadians = 3.9269908169872414;
      frDetector->yArmAltitudeRadians = 0.;
      frDetector->yArmAzimuthRadians = 5.497787143782138;
      detector = XLALCreateDetector(network + k, frDetector, LALDETECTORTYPE_IFODIFF);
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


/* Compute horizon distance for a given detector.  After equation (2) from LLOID paper. */
static double horizon_distance(double m1, double m2, double f_low, double f_high, double snr_thresh, char detector)
{
  double (*S)(double);
  double Dhor, f;

  /* Chirp mass in seconds */
  const double Mchirp = LAL_MTSUN_SI * pow(pow(m1 * m2, 3) / (m1 + m2), 0.2);

  /* Integration step in Hz */
  static const double df = 0.01;

  /* Select the function to use for sampling the PSD */
  switch (detector)
    {
    case 'H':
    case 'L':
    case 'I':
      S = XLALSimNoisePSDaLIGOZeroDetHighPower;
    break;
    case 'V':
      S = XLALSimNoisePSDAdvVirgo;
      break;
    case 'K':
      S = XLALSimNoisePSDKAGRA;
      break;
    default:
      fprintf(stderr, "This line should never be reached.\n");
      exit(EXIT_FAILURE);
    }

  for (Dhor = 0, f = f_low; f < f_high; f += df)
    Dhor += pow(f, -7./3) / S(f) * df;

  Dhor = sqrt(5./6 * Dhor);
  Dhor *= pow(Mchirp, 5./6) * pow(M_PI, -2./3) * LAL_C_SI / snr_thresh;

  /* Convert from m to Mpc */
  Dhor /= 1e6 * LAL_PC_SI;

  return Dhor;
}
