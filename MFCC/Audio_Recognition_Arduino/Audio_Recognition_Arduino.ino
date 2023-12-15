#include <Arduino.h>
#include <Serial.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

// Import TensorFlow stuff
#include "TensorFlowLite.h"
#include "tensorflow/lite/micro/all_ops_resolver.h"
#include "tensorflow/lite/micro/micro_error_reporter.h"
#include "tensorflow/lite/micro/micro_interpreter.h"
#include "tensorflow/lite/version.h"

// Our model
#include "model_tflite.h"

// Other headers
#include "audio_data.h"
#include "fftmat.h"
#include "window.h"

#define WINDOW_SIZE 64
#define FFTPOINT 64
#define HOPSIZE 32
#define TOTAL_SAMPLES 48000
#define NUM_FILTERS 40
#define SAMPLE_RATE 16000

#define OUTLEN 13
#define BAUD_RATE 9600
#define LED 22

float xn[TOTAL_SAMPLES];
float y[OUTLEN];

const int TLEN = (TOTAL_SAMPLES - WINDOW_SIZE) / HOPSIZE + 1;
const int FLEN = FFTPOINT / 2 + 1;

float stft[TLEN * FLEN];
float x_mel[TLEN * NUM_FILTERS];

float x_window[FFTPOINT];
float xfft_r[FFTPOINT];
float xfft_i[FFTPOINT];
float mult_window[FFTPOINT];

void dftcustom( int N);
int stftcustom( int window_size , int hopsize , int fftsize);
void pre_emphasis(float* inp, int len, float alpha);
void power_spectrum(float* inp, int time_len, int fft_len);
void mfcc(float* inp, float* out, int inp_len, int fft_len, int num_filters, float sr);


// Figure out what's going on in our model
#define DEBUG 0

// TFLite globals, used for compatibility with Arduino-style sketches
namespace {
  tflite::ErrorReporter* error_reporter = nullptr;
  const tflite::Model* model = nullptr;
  tflite::MicroInterpreter* interpreter = nullptr;
  TfLiteTensor* model_input = nullptr;
  TfLiteTensor* model_output = nullptr;

  // Create an area of memory to use for input, output, and other TensorFlow
  // arrays. You'll need to adjust this by combiling, running, and looking
  // for errors.
  constexpr int kTensorArenaSize = 50 * 1024;
  uint8_t tensor_arena[kTensorArenaSize];
} // namespace


void setup() {

  pinMode(LED, OUTPUT);
  digitalWrite(LED, HIGH);

  Serial.begin(BAUD_RATE);  
  while (!Serial);

  // Set up logging (will report to Serial, even within TFLite functions)
  static tflite::MicroErrorReporter micro_error_reporter;
  error_reporter = &micro_error_reporter;

  // Map the model into a usable data structure
  model = tflite::GetModel(model_tflite);
  if (model->version() != TFLITE_SCHEMA_VERSION) {
    error_reporter->Report("Model version does not match Schema");
    while(1);
  }

  // Pull in only needed operations (should match NN layers)
  // Available ops:
  //  https://github.com/tensorflow/tensorflow/blob/master/tensorflow/lite/micro/kernels/micro_ops.h
  
  static tflite::AllOpsResolver resolver;

  // Build an interpreter to run the model
  static tflite::MicroInterpreter static_interpreter(
    model, resolver, tensor_arena, kTensorArenaSize, error_reporter);
  interpreter = &static_interpreter;

  // Allocate memory from the tensor_arena for the model's tensors
  TfLiteStatus allocate_status = interpreter->AllocateTensors();
  if (allocate_status != kTfLiteOk) {
    error_reporter->Report("AllocateTensors() failed");
    while(1);
  }

  // Assign model input and output buffers (tensors) to pointers
  model_input = interpreter->input(0);
  model_output = interpreter->output(0);

 delay(500);
  if(DEBUG) {
    Serial.print("Number of input dimensions: ");
    Serial.println(model_input->dims->size);
    Serial.print("Dim 1 size: ");
    Serial.println(model_input->dims->data[0]);
    Serial.print("Dim 2 size: ");
    Serial.println(model_input->dims->data[1]);
    Serial.print("Dim 3 size: ");
    Serial.println(model_input->dims->data[2]);
    Serial.print("Dim 4 size: ");
    Serial.println(model_input->dims->data[3]);
    Serial.print("Input type: ");
    Serial.println(model_input->type);
    Serial.print("Number of output dimensions: ");
    Serial.println(model_output->dims->size);
    Serial.print("Dim 1 size: ");
    Serial.println(model_output->dims->data[0]);
    Serial.print("Dim 2 size: ");
    Serial.println(model_output->dims->data[1]);
    Serial.print("Output type: ");
    Serial.println(model_output->type);
  }
  delay(500);
}

void loop() {
  // put your main code here, to run repeatedly

  for(int i=0; i < TOTAL_SAMPLES; i++)
    xn[i] = 0;
  
  while(Serial.read() != 'S');          // send request
  Serial.write('Y');
  
  digitalWrite(LED, LOW);

  while(Serial.read() != 'f');
  int amplifier = Serial.parseInt();

 

  // Performe MFCC
  pre_emphasis(audioData, TOTAL_SAMPLES, 0.97);
  stftcustom(WINDOW_SIZE, HOPSIZE, FFTPOINT);
  power_spectrum(stft, TLEN * FLEN, FFTPOINT);
  mfcc(stft, x_mel, TLEN, FFTPOINT, NUM_FILTERS, SAMPLE_RATE);

  // Model input
  for(int i = 0; i < TLEN * NUM_FILTERS; i++){
    model_input->data.f[i] = x_mel[i];
  }

  // Start Timer
  unsigned long start_time = millis();
  
  // Invoke
  TfLiteStatus invoke_status = interpreter->Invoke();
  if (invoke_status != kTfLiteOk) {
    error_reporter->Report("Invoke failed!");
  }

  // End Timer
  unsigned long end_time = millis();
  unsigned long elapsed_time = end_time - start_time;

  // Print MFCC coefficients directly to the terminal
for(int i = 0; i < OUTLEN; i++)
{
    y[i] = model_output->data.int8[i];
    Serial.print(y[i]);
    Serial.print(" ");
}

// Print a newline character to separate the next set of MFCC coefficients
Serial.println();

  Serial.println(elapsed_time);

  digitalWrite(LED, HIGH);

}


// -----------------------------------------------------------------------------------------------------------------------------------------------
/* Define required functions here */


#define FFT_LEN 64

// Slightly amplify higher frequency signals to compensate for human speech characteristics (low freq has higher energy)
void pre_emphasis(float* inp, int len, float alpha) {
    float prev_inp = 0;
    float temp;

    for(int i = 0; i < len; i++) {
        temp        = inp[i];                       // Save current input for later
        inp[i]      = inp[i] - alpha * prev_inp;    // In-place FIR filtering
        prev_inp    = temp;                         // Current input for next sample
    }
}

void dftcustom( int N) {
   int i=0, j;
   float temp;
   float a, b;

   for (i = 0; i < (N/2 + 1); i++){      
      temp = 0;
      for (j = 0; j < N; j++)
        temp = temp + mult_window[j]*Rmat[i][j];
      xfft_r[i] = temp;

      temp = 0;
      for (j = 0;j < N; j++)
        temp = temp + mult_window[j]*Imat[i][j];
      xfft_i[i] = temp;
   }
}

int stftcustom( int window_size , int hopsize , int fftsize )
{
    int start_index = 0-hopsize;
    int end_index = window_size - hopsize;
    int counter=0;
    int total_samples =  TOTAL_SAMPLES;
    int length;
    int i, j;
    int timebins=0;

    while((start_index + hopsize) <= (total_samples - window_size))
    {
        memset(x_window, 0, sizeof(x_window));
        memset(mult_window, 0, sizeof(mult_window));

        start_index = start_index + hopsize;
        end_index = end_index + hopsize;

        for (i = start_index; i < end_index; i++)
            x_window[i-start_index] = audioData[i];

        length = end_index - start_index;

        for (i=0;i<length;i++)
            mult_window[i] =  x_window[i] * window_fn[i];

        dftcustom(length);

        for(i = 0; i < FLEN; i++)
            stft[timebins * FLEN + i] = sqrt(xfft_r[i]*xfft_r[i] + xfft_i[i]*xfft_i[i]);

        timebins++;
    }
    return timebins;
}

// Convert STFT output to power domain by squaring and scaling
void power_spectrum(float* inp, int inp_len, int fft_len) {
    float fac = 1.0 / (float)fft_len;

    for(int i = 0; i < inp_len; i++)
        inp[i] = pow(inp[i], 2) * fac;
}

// Helper function for mel filtering
static inline float get_mel_from_hertz(float hertz) {
    return 2595 * log10(1 + (hertz/ 700));
}

// Helper function for mel filtering
static inline float get_hertz_from_mel(float mel) {
    return 700 * (pow(10, mel / 2595) - 1);
}

// Perform mfcc on whole data
void mfcc(float* inp, float* out, int inp_len, int fft_len, int num_filters, float sr) {
    int inp_width = 1 + fft_len / 2;
    float max_mel = get_mel_from_hertz(sr / 2);
    float sum;

    for(int i = 0; i < num_filters; i++) {
      
      float start_hz = get_hertz_from_mel(max_mel * (((float)  i  ) / ((float)num_filters + 2)));
      float peak_hz  = get_hertz_from_mel(max_mel * (((float)i + 1) / ((float)num_filters + 2)));
      float end_hz   = get_hertz_from_mel(max_mel * (((float)i + 2) / ((float)num_filters + 2)));
      
      int start_idx  = floor(start_hz * (fft_len + 1) / sr);
      int peak_idx   = floor(peak_hz * (fft_len + 1) / sr);
      int end_idx    = floor(end_hz * (fft_len + 1) / sr);
      
      for(int j = 0; j < inp_len; j++) {
        out[i * inp_len + j] = inp[j * inp_width + peak_idx];
        
        // Ascending
        for(int k = start_idx; k < peak_idx; k++) 
            out[i * inp_len + j] += inp[j * inp_width + k] * (k - start_idx) / (peak_idx - start_idx);

        // Descending
        for(int k = peak_idx + 1; k < end_idx; k++) 
            out[i * inp_len + j] += inp[j * inp_width + k] * (end_idx - k) / (end_idx - peak_idx);

      }
   }
}
