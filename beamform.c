// 3D Ultrasound beamforming baseline code for EECS 570 
// Created by: Richard Sampson, Amlan Nayak, Thomas F. Wenisch
// Revision 1.0 - 11/15/16

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include "sys/time.h"
#include <pthread.h>

#define NUM_THREADS 10
#define NUM_THREADS2 1024
pthread_mutex_t lock=PTHREAD_MUTEX_INITIALIZER;

struct thread_args{
    int start;
    int end;
};

struct thread_args2{
    int start;
    int end;
    int rx;
};

int boundry;
int outer_boundry;
float *point_x; // Point x position
float *point_y; // Point y position
float *point_z; // Point z position
float tx_x = 0; // Transmit transducer x position
float tx_y = 0; // Transmit transducer y position
float tx_z = -0.001; // Transmit transducer z position
float *dist_tx; // Transmit distance (ie first leg only)
int it_rx; // Iterator for recieve transducer
float *image;  // Pointer to full image (accumulated so far)
const float idx_const = 0.000009625; // Speed of sound and sampling rate, converts dist to index
const int filter_delay = 140; // Constant added to index to account filter delay (off by 1 from MATLAB)
int data_len = 12308; // Number for pre-processed data values per channel
float *rx_x; // Receive transducer x position
float *rx_y; // Receive transducer y position
float rx_z = 0; // Receive transducer z position
float *rx_data; // Pointer to pre-processed receive channel data

void *cal_tx_data(void *args){
	struct thread_args *range = (struct thread_args *) args;
	float x_comp; // Itermediate value for dist calc
	float y_comp; // Itermediate value for dist calc
	float z_comp; // Itermediate value for dist calc
	int point = 0;
	for(point = range->start; point < range->end; point++){

		x_comp = tx_x - point_x[point];
		x_comp = x_comp * x_comp;
		y_comp = tx_y - point_y[point];
		y_comp = y_comp * y_comp;
		z_comp = tx_z - point_z[point];
		z_comp = z_comp * z_comp;
		dist_tx[point] = (float)sqrt(x_comp + y_comp + z_comp);
	}
	pthread_exit(0);
}

void *cal_rx_data(void *args){

		float *image_pos = image; // Reset image pointer back to beginning
		int offset = 0;
		int point = 0; // Reset 
		int ppoint;
		int bboundry = 64*64*1560;
		struct thread_args2 *range = (struct thread_args2 *) args;
		float x_comp; // Itermediate value for dist calc
		float y_comp; // Itermediate value for dist calc
		float z_comp; // Itermediate value for dist calc
		float dist; // Full distance
		int index; // Index into transducer data
		int rx=0;
		// Iterate over entire image space
		for (point = range->start; point < range->end; point++ ){
			if((point+range->rx)>=6389760){
				ppoint = (point+range->rx)-6389760; 
			}
			else {
				ppoint = point  + range->rx;
								
			}
			
			rx = range->rx;
			offset = rx *data_len;

			x_comp = rx_x[rx] - point_x[ppoint];
			x_comp = x_comp * x_comp;
			y_comp = rx_y[rx] - point_y[ppoint];
			y_comp = y_comp * y_comp;
			z_comp = rx_z - point_z[ppoint];
			z_comp = z_comp * z_comp;

			long thread_private_tmp = 0;
			dist = dist_tx[ppoint] + (float)sqrt(x_comp + y_comp + z_comp);
			index = (int)(dist/idx_const + filter_delay + 0.5);
            image_pos[ppoint]  += rx_data[index+offset];

		}
	pthread_exit(0);
}


int main (int argc, char **argv) {

	int size = atoi(argv[1]);

	/* Variables for transducer geometry */
	int trans_x = 32; // Transducers in x dim
	int trans_y = 32; // Transducers in y dim
	
	int offset = 0; // Offset into rx_data



	/* Variables for image space points */
	int point; // Index into image space


	int pts_r = 1560; // Radial points along scanline
	int sls_t = size; // Number of scanlines in theta
	int sls_p = size; // Number of scanlines in phi

	float *image_pos; // Pointer to current position in image

	/* Iterators */
	int it_r; // Iterator for r
	int it_t; // Iterator for theta
	int it_p; // Iterator for phi

	/* Variables for distance calculation and index conversion */
	float x_comp; // Itermediate value for dist calc
	float y_comp; // Itermediate value for dist calc
	float z_comp; // Itermediate value for dist calc

	//float *dist_tx; // Transmit distance (ie first leg only)

        FILE* input;
        FILE* output;

	/* Allocate space for data */
	rx_x = (float*) malloc(trans_x * trans_y * sizeof(float));
	if (rx_x == NULL) fprintf(stderr, "Bad malloc on rx_x\n");
	rx_y = (float*) malloc(trans_x * trans_y * sizeof(float));
	if (rx_y == NULL) fprintf(stderr, "Bad malloc on rx_y\n");
	rx_data = (float*) malloc(data_len * trans_x * trans_y * sizeof(float));
	if (rx_data == NULL) fprintf(stderr, "Bad malloc on rx_data\n");

	point_x = (float *) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (point_x == NULL) fprintf(stderr, "Bad malloc on point_x\n");
	point_y = (float *) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (point_y == NULL) fprintf(stderr, "Bad malloc on point_y\n");
	point_z = (float *) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (point_z == NULL) fprintf(stderr, "Bad malloc on point_z\n");

	dist_tx = (float*) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (dist_tx == NULL) fprintf(stderr, "Bad malloc on dist_tx\n");

	image = (float *) malloc(pts_r * sls_t * sls_p * sizeof(float));
	if (image == NULL) fprintf(stderr, "Bad malloc on image\n");
	memset(image, 0, pts_r * sls_t * sls_p * sizeof(float));

	/* validate command line parameter */
	if (argc < 1 || !(strcmp(argv[1],"16") || strcmp(argv[1],"32") || strcmp(argv[1],"64"))) {
	  printf("Usage: %s {16|32|64}\n",argv[0]);
	  fflush(stdout);
	  exit(-1);
	}

	char buff[128];
        #ifdef __MIC__
	  sprintf(buff, "./beamforming_input_%s.bin", argv[1]);
        #else // !__MIC__
	  sprintf(buff, "./beamforming_input_%s.bin", argv[1]);
	  //sprintf(buff, "/n/typhon/data1/home/eecs570/beamforming_input_%s.bin", argv[1]);
        #endif

        input = fopen(buff,"rb");
	if (!input) {
	  printf("Unable to open input file %s.\n", buff);
	  fflush(stdout);
	  exit(-1);
	}	

	/* Load data from binary */
	fread(rx_x, sizeof(float), trans_x * trans_y, input); 
	fread(rx_y, sizeof(float), trans_x * trans_y, input); 

	fread(point_x, sizeof(float), pts_r * sls_t * sls_p, input); 
	fread(point_y, sizeof(float), pts_r * sls_t * sls_p, input); 
	fread(point_z, sizeof(float), pts_r * sls_t * sls_p, input); 

	fread(rx_data, sizeof(float), data_len * trans_x * trans_y, input); 
        fclose(input);

	printf("Beginning computation\n");
	fflush(stdout);

	/* get start timestamp */
 	struct timeval tv;
    	gettimeofday(&tv,NULL);
    	uint64_t start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
 
	/* --------------------------- COMPUTATION ------------------------------ */
	/* First compute transmit distance */
	boundry 	= sls_t*sls_p*pts_r;

	pthread_t 		child_threads1[NUM_THREADS];
    	struct thread_args 	work_ranges1[NUM_THREADS];
    	long int 			current_start, range;
    	long int range2;
	int i,j,k,l,m,n;
    	current_start = 0;
    	range = boundry / NUM_THREADS;
    	for(i = 0; i < NUM_THREADS; i++) {
    	    work_ranges1[i].start = current_start;
    	    work_ranges1[i].end = current_start + range;
    	    current_start += range;
    	}
    	work_ranges1[NUM_THREADS-1].end = boundry;

	for(i = 0; i < NUM_THREADS; i++) {
        	pthread_create(&child_threads1[i], NULL, cal_tx_data, &work_ranges1[i]);
    	}
    	for(i = 0; i < NUM_THREADS; i++) {
        	pthread_join(child_threads1[i], NULL);
    	}

	/* Now compute reflected distance, find index values, add to image */
		pthread_t 		child_threads2[NUM_THREADS2];
    	struct thread_args2 	work_ranges2[NUM_THREADS2];
    	current_start = 0;
    	range2 =  boundry / NUM_THREADS2;
    	range2 =  range2 * trans_x*trans_y;
    	int rx=0;
       	for(j = 0; j < NUM_THREADS2; j++) {
    	    work_ranges2[j].start = current_start;
    	    work_ranges2[j].end = current_start + range2;
    	    work_ranges2[j].rx = rx; 
			rx += 1 ;
    	    current_start =	work_ranges2[j].end % boundry;
    	    printf("start is %d\n",work_ranges2[j].rx);
    	    printf("start is %d\n",work_ranges2[j].start);
    	    printf("start is %d\n",work_ranges2[j].end);
    	    printf("===============================\n");
    	}

    	work_ranges2[NUM_THREADS2-1].end = range2;

		for(i = 0; i < NUM_THREADS2; i++) {
        		pthread_create(&child_threads2[i], NULL, cal_rx_data, &work_ranges2[i]);
    		}
    		for(i = 0; i < NUM_THREADS2; i++) {
        		pthread_join(child_threads2[i], NULL);
    		}

	

	/* --------------------------------------------------------------------- */

	/* get elapsed time */
    	gettimeofday(&tv,NULL);
    	uint64_t end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    	uint64_t elapsed = end - start;

	printf("@@@ Elapsed time (usec): %lld\n", elapsed);
	printf("Processing complete.  Preparing output.\n");
	fflush(stdout);

	/* Write result to file */
	char* out_filename;
        #ifdef __MIC__
	  out_filename = "/home/micuser/beamforming_output.bin";
        #else // !__MIC__
	  out_filename = "beamforming_output.bin";
        #endif
        output = fopen(out_filename,"wb");
	fwrite(image, sizeof(float), pts_r * sls_t * sls_p, output); 
	fclose(output);

	printf("Output complete.\n");
	fflush(stdout);

	/* Cleanup */
	free(rx_x);
	free(rx_y);
	free(rx_data);
	free(point_x);
	free(point_y);
	free(point_z);
	free(dist_tx);
	free(image);

	return 0;
}
