// 3D Ultrasound beamforming baseline code for EECS 570 
// Created by: Richard Sampson, Amlan Nayak, Thomas F. Wenisch
// Revision 1.0 - 11/15/16

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <pthread.h>

#define NUM_THREADS 300

pthread_mutex_t lock=PTHREAD_MUTEX_INITIALIZER;

struct thread_args{
    int start;
    int end;
};


int boundry;
float *point_x; // Point x position
float *point_y; // Point y position
float *point_z; // Point z position
float tx_x = 0; // Transmit transducer x position
float tx_y = 0; // Transmit transducer y position
float tx_z = -0.001; // Transmit transducer z position
float *dist_tx; // Transmit distance (ie first leg only)
float *image;  // Pointer to full image (accumulated so far)
const float idx_const = 0.000009625; // Speed of sound and sampling rate, converts dist to index
const int filter_delay = 140; // Constant added to index to account filter delay (off by 1 from MATLAB)
int data_len = 12308; // Number for pre-processed data values per channel
float *rx_x; // Receive transducer x position
float *rx_y; // Receive transducer y position
float rx_z = 0; // Receive transducer z position
float *rx_data; // Pointer to pre-processed receive channel data

int trans_x = 32; // Transducers in x dim
int trans_y = 32; // Transducers in y dim

void *cal_dist(void *args){
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

void *cal_image(void *args){
		struct thread_args *range = (struct thread_args *) args;
		float *image_pos = image; // Reset image pointer back to beginning
		int offset = 0;
		int point = 0; // Reset    		
		float x_comp; // Itermediate value for dist calc
		float y_comp; // Itermediate value for dist calc
		float z_comp; // Itermediate value for dist calc
		float dist; // Full distance
		int index; // Index into transducer data
		int j;
		
		for (j=0;j <trans_x*trans_y ;j++){
	
		for (point = range->start; point < range->end; point++ ){	

			x_comp = rx_x[j] - point_x[point];
			x_comp = x_comp * x_comp;
			y_comp = rx_y[j] - point_y[point];
			y_comp = y_comp * y_comp;
			z_comp = rx_z - point_z[point];
			z_comp = z_comp * z_comp;

			dist = dist_tx[point] + (float)sqrt(x_comp + y_comp + z_comp);
			index = (int)(dist/idx_const + filter_delay + 0.5);
            		image_pos[point]  += rx_data[index+offset];

		}
		offset += data_len;
		}
	pthread_exit(0);
}


int main (int argc, char **argv) {

	int size = atoi(argv[1]);

	/* Variables for transducer geometry */


	int pts_r = 1560; // Radial points along scanline
	int sls_t = size; // Number of scanlines in theta
	int sls_p = size; // Number of scanlines in phi


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
	  sprintf(buff, "/beamforming_input_%s.bin", argv[1]);
        #else // !__MIC__
	  sprintf(buff, "/n/typhon/data1/home/eecs570/beamforming_input_%s.bin", argv[1]);
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

	pthread_t 		child_threads[NUM_THREADS];
    	struct thread_args 	work_ranges[NUM_THREADS];
    	long int 			current_start, range;
	int i;
    	current_start = 0;
    	range = boundry / NUM_THREADS;
    	for(i = 0; i < NUM_THREADS; i++) {
    	    work_ranges[i].start = current_start;
    	    work_ranges[i].end = current_start + range;
    	    current_start += range;
    	}
    	work_ranges[NUM_THREADS-1].end = boundry;

	for(i = 0; i < NUM_THREADS; i++) {
        	pthread_create(&child_threads[i], NULL, cal_dist, &work_ranges[i]);
    	}
    	for(i = 0; i < NUM_THREADS; i++) {
        	pthread_join(child_threads[i], NULL);
    	}


	for(i = 0; i < NUM_THREADS; i++) {
		pthread_create(&child_threads[i], NULL, cal_image, &work_ranges[i]);
    	}
    	for(i = 0; i < NUM_THREADS; i++) {
        	pthread_join(child_threads[i], NULL);
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
