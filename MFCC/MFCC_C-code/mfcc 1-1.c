#include <stdio.h>
#include "stdlib.h"
#include <math.h>

#define DEBUG 1

#define WINDOW_SIZE 500
#define HOP_SIZE (WINDOW_SIZE/4)
#define FFT_SIZE 512
#define N FFT_SIZE
#define TOTAL_SAMPLES 104429  // we can declare this in the audio header file
#define STFT_COL  ((FFT_SIZE/2) + 1)
#define STFT_ROWS  (((TOTAL_SAMPLES - WINDOW_SIZE)/HOP_SIZE) + 1)

#include "audio_sample.h"   //only to be used during debug
#include "fft_table.h"
//#include<mel_freq_bins.h> //to be included which would contial the nearest FFT bins
//#include <dct_cosine_table.h>                //This LUT needs to be include for DCT  

//MFCC parameters
#define ALPHA 0.97
#define SAMPLE_RATE 22050
#define NO_OF_FILTERS 40
#define MACHINE_EPSILON (2.22044604925e-16)    //Machine_Epsilon
#define PI 3.14159265359

void dft(float samples[], float dft_matrix_real[], float dft_matrix_imag[])
{
    int twf_exp;
    float sum_real = 0, sum_imag = 0;
    for (int k=0; k<((FFT_SIZE/2) +1); k++)
    {
        for(int n=0; n<FFT_SIZE; n++)
        {
            twf_exp = (n*k) % FFT_SIZE;
            if(twf_exp < FFT_SIZE/2)
            {
                sum_real = sum_real + (cos_table[twf_exp] * samples[n]) ;
                sum_imag = sum_imag + (sin_table[twf_exp] * samples[n]) ;
            }
            else
            {
                twf_exp = twf_exp - (FFT_SIZE/2);
                sum_real = sum_real + (cos_table[twf_exp] * samples[n] * (-1));
                sum_imag = sum_imag + (sin_table[twf_exp] * samples[n] * (-1));
            }
            // printf("sample : %f\n", samples[n]);
        }
        dft_matrix_real[k] = sum_real;
        dft_matrix_imag[k] = sum_imag;
        sum_real = 0; sum_imag = 0;
    }

}
    


void stft(float stft_result_real[STFT_ROWS][STFT_COL], float stft_result_imag[STFT_ROWS][STFT_COL])
{
    int total_samples = TOTAL_SAMPLES ;
    int start_index = 0 - HOP_SIZE;    //int end_index = WINDOW_SIZE - HOP_SIZE;
    // float *stft_real = NULL, *stft_imag = NULL;
    int counter = 0, diff = FFT_SIZE - WINDOW_SIZE;
    float x_window[FFT_SIZE], mult_window[FFT_SIZE];
    float dft_matrix_real[(FFT_SIZE/2) + 1];
    float dft_matrix_imag[(FFT_SIZE/2) + 1];
    while(start_index + HOP_SIZE + WINDOW_SIZE <= total_samples)
    {
        start_index = start_index + HOP_SIZE;
        //end_index = end_index + HOP_SIZE;

        //Making the window size equal to FFT size, by paddind zeroes on both end
        for(int i=0; i<(diff/2); i++)
            x_window[i] = 0;
        for(int i=0; i<(WINDOW_SIZE); i++)
        {
            x_window[(diff/2) + i] = audio_samples[start_index+i];
            // printf("\n x_window_before : %e", x_window[i]);
        }
        for(int i = ((diff/2)+WINDOW_SIZE); i<FFT_SIZE ; i++)
            x_window[i] = 0;

        //Multiplying window_function with x_window;
        for(int i=0; i<FFT_SIZE; i++)
        {
            x_window[i] = x_window[i] * window_func[i];
            // printf("\n x_window_before : %e", x_window[i]);

        }
        
        //do STFT on the window function
        dft(x_window, dft_matrix_real, dft_matrix_imag);
        for(int i=0; i<STFT_COL; i++)
        {
            stft_result_real[counter][i] = dft_matrix_real[i];
            stft_result_imag[counter][i] = dft_matrix_imag[i];
        }
        counter = counter + 1;
        //if(counter==700)
        //    for(int i =0; i<(FFT_SIZE/2 + 1); i++)
        //        printf("real : %E, imag : %e\n", dft_matrix_real[i], dft_matrix_imag[i]);
    }
    #ifdef DEBUG
        printf("End of STFT\n");
    #endif
    
}

void pre_emphasis()
{   
    #ifdef DEBUG
        printf("In pre_emphasis\n");
    #endif
    for(int i=TOTAL_SAMPLES-1; i>=1; i--)
    {
        audio_samples[i] = audio_samples[i] - (ALPHA * (audio_samples[i-1]));
    }
}

float get_mel_from_hz(float hertz)
{
    return (float)(2595 * (log10(1 + (hertz/700))));
}

float get_hz_from_mel(float mel)
{
    return(float)(700 * ((pow(10, (mel/2595))) - 1));
}

void get_triangular_filters(float prev_bin, float center_bin, float next_bin, int bin_fb, float filter_banks[NO_OF_FILTERS][STFT_COL])
{
    //positive slope side
    for(int freq=prev_bin; freq<center_bin; freq++)
    {
        filter_banks[bin_fb-1][freq] = (freq - prev_bin)/(center_bin - prev_bin);
    }
    //Negtive slope
    for(int freq = center_bin+1; freq<next_bin; freq++)
    {
        filter_banks[bin_fb-1][freq] = (next_bin - freq)/(next_bin - center_bin);
    }
    //Tip of the filter
    filter_banks[bin_fb-1][(int)center_bin] = 1;
    //Replacing zeroes in the mel_filters with machine_epsilon
    for(int j =0; j<STFT_COL; j++)
    {
        if(filter_banks[bin_fb-1][j]==0)
            filter_banks[bin_fb-1][j] == MACHINE_EPSILON;
    }
}

void mel_filter_bank(float stft_result_real[STFT_ROWS][STFT_COL], float filter_banks[NO_OF_FILTERS][STFT_COL])
{
    float min_mel =0;
    float max_mel = get_mel_from_hz(SAMPLE_RATE/2);
    float mel_filter_points[NO_OF_FILTERS + 2];
    float mel_filter_spacing = max_mel/(NO_OF_FILTERS + 1);
    //computing mel frequencies
    for(int i=0; i<(NO_OF_FILTERS +2); i++)
    {
        if(i==0)
            mel_filter_points[i] = 0;
        else
        {
            mel_filter_points[i] = mel_filter_points[i-1] + mel_filter_spacing;
        }
        #ifdef DEBUG
            printf("Mel filter point %d: %f", i, mel_filter_points[i]);
            printf("\n");
        #endif
    }
    //Converting mel frequencies to hz and storing it in the variable 
    for(int i=1; i<(NO_OF_FILTERS + 2); i++)
    {
        mel_filter_points[i] = get_hz_from_mel(mel_filter_points[i]);
        #ifdef DEBUG
            printf("Mel to hz: %f \n", mel_filter_points[i]);
        #endif
    }
    //connverting the frequencies to the nearest bin
    for(int i=0; i<(NO_OF_FILTERS + 2); i++)
    {
        mel_filter_points[i] = (((int)mel_filter_points[i]) * (FFT_SIZE + 1))/(SAMPLE_RATE);
        #ifdef DEBUG
            printf("FFT_bins : %f \n", mel_filter_points[i]);
        #endif
    }
    //genarating triangular filters
    for(int bin_fb=1; bin_fb<(NO_OF_FILTERS + 1); bin_fb++)
    {
        float prev_bin = mel_filter_points[bin_fb-1];
        float center_bin = mel_filter_points[bin_fb];
        float next_bin = mel_filter_points[bin_fb+1];
        get_triangular_filters(prev_bin, center_bin, next_bin, bin_fb, filter_banks);
       
        #ifdef DEBUG
        printf("Triangulat filter %d values:", bin_fb-1);
        for(int i=prev_bin; i<=next_bin; i++)
            printf("  %e", filter_banks[bin_fb-1][i]);
        printf("\n");
        #endif
    }

}

int main()
{   
    //this two arrays store the stft result
    float stft_result_real[STFT_ROWS][STFT_COL], stft_result_imag[STFT_ROWS][STFT_COL];

    //to store filter banks
    float filter_banks[NO_OF_FILTERS][STFT_COL];

    #ifdef DEBUG
        //for(int i=1000;i<2000;i++)
            printf("1999 :%e \n",audio_samples[1999]);
            printf("2000: %e \n",audio_samples[2000]);
    #endif


    // this function will provide the pre_emphasied audio signal
    pre_emphasis();
    #ifdef DEBUG
        //for(int i=1000;i<2000;i++)
            printf("2000: %e \n",audio_samples[2000]);
    #endif

    // this will calculate the STFT of the preemphasized audio signal
    stft(stft_result_real, stft_result_imag);
    #ifdef DEBUG
        for(int i =0; i<(STFT_COL); i++)
        {
            printf("\nSTFT_real : %e, imag : %e\n", stft_result_real[770][i], stft_result_imag[770][i]);
        }
    #endif
    // //Calculating power at each frequency, i.e geneting power spectrum
    // // the power spectrum is stored in the variable stft_result_real[STFT_ROWS][STFT_COL]
    #ifdef DEBUG
        printf("\nPower Spectrum\n");
    #endif
    for(int i=0; i<STFT_ROWS; i++)
    {
        for(int j=0; j<STFT_COL; j++)
        {
            stft_result_real[i][j] = ((stft_result_real[i][j] * stft_result_real[i][j]) + (stft_result_imag[i][j] * stft_result_imag[i][j])) / FFT_SIZE ;
             #ifdef DEBUG
                if(i==770)
                    printf("%e \n", stft_result_real[i][j]);
            #endif
        }
    }

    //creating mel center frequencies

    for(int i=0; i <NO_OF_FILTERS; i++)
    {
        for(int j=0; j<STFT_COL; j++)
            filter_banks[i][j] = 0;
    }
    mel_filter_bank(stft_result_real, filter_banks);
    
    //Taking dot product of mel_filters and power spectrograms
    //storing result in stft_result_imag
    for(int window = 0; window<STFT_ROWS; window++)
    { 
        for(int filter=0; filter<NO_OF_FILTERS; filter++)
        {
            float sum = 0;
            for(int i=0; i<STFT_COL; i++)
            {
                sum = sum + (filter_banks[filter][i] * stft_result_real[window][i]); 
            }
            stft_result_imag[window][filter] = sum;
            //Adding machine_epsilon to each element and takin log
            stft_result_imag[window][filter] = stft_result_imag[window][filter] + MACHINE_EPSILON;
            stft_result_imag[window][filter] = (float)log(stft_result_imag[window][filter]);
        }
    }
    
    //Now the dot product of mel filters and spectogram is stored in
    //stft_result_real[STFT_ROWS][NO_OF_FILTERS]

    #ifdef DEBUG
    printf("\nAfter applying Mel filter banks on spectrogram and taking natural log:\n ");
        for(int i=0; i<NO_OF_FILTERS; i++)
        {
            printf("  %f", stft_result_imag[770][i]);
        }
    printf("\n");
    #endif

    //Doing DCT
    for(int window=0; window<STFT_ROWS; window++)
    {
       
        for(int k=0; k<(NO_OF_FILTERS); k++)
        {
            float sum=0;
            for(int n=0; n<(NO_OF_FILTERS); n++)
            {   
                float cos_numerator = (2 * PI * k *n) + (PI * k);
                float cos_denomenator = 2 * (NO_OF_FILTERS);
                sum = sum + (stft_result_imag[window][n] * ((float)cos(cos_numerator/cos_denomenator))); 
            }
            stft_result_real[window][k]= 2*sum;
        }
    }
    #ifdef DEBUG
    printf("\nAfter applying DCT\n ");
    for(int i=0; i<NO_OF_FILTERS; i++)
    {
        printf("  %f", stft_result_real[770][i]);
    }
    printf("\n");
    #endif

    
    return 0;
}