///my definitions
#define TOTALFRAMES NUM_SECONDS*SAMPLE_RATE
#define NUMSAMPLES TOTALFRAMES*NUM_CHANNELS
#define TEST_PRINTS 1
#define IS_UBUNTU 1
#define USE_FRAME 1
#define USE_DTW 1

#if USE_DTW
#include "dtw.h"
#define DTW_PASS_VALUE 4
#define COEFFFS 12
#else
#define COEFFFS 100
#endif

#if IS_UBUNTU
#include <termios.h>
#include <unistd.h>
/* reads from keypress, doesn't echo */
int getch(void)
{
    struct termios oldattr, newattr;
    int ch;
    tcgetattr( STDIN_FILENO, &oldattr );
    newattr = oldattr;
    newattr.c_lflag &= ~( ICANON | ECHO );
    tcsetattr( STDIN_FILENO, TCSANOW, &newattr );
    ch = getchar();
    tcsetattr( STDIN_FILENO, TCSANOW, &oldattr );
    return ch;
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "libmfcc.h"


///portaudio definitions
 /* #define SAMPLE_RATE (44100) // Test failure to open with this value. */
 #define SAMPLE_RATE (44100)
 #define FRAMES_PER_BUFFER (512)
 #define NUM_SECONDS (2)
 #define NUM_CHANNELS (2)
 /* #define DITHER_FLAG (paDitherOff) */
 #define DITHER_FLAG (0)

 #define WRITE_TO_FILE (0)

 /* Select sample format. */
 #if 1
 #define PA_SAMPLE_TYPE paFloat32
 typedef float SAMPLE;
 #define SAMPLE_SILENCE (0.0f)
 #define PRINTF_S_FORMAT "%.8f"
 #elif 1
#define PA_SAMPLE_TYPE paInt16
 typedef short SAMPLE;
 #define SAMPLE_SILENCE (0)
#define PRINTF_S_FORMAT "%d"
 #elif 0
 #define PA_SAMPLE_TYPE paInt8
 typedef char SAMPLE;
 #define SAMPLE_SILENCE (0)
 #define PRINTF_S_FORMAT "%d"
 #else
 #define PA_SAMPLE_TYPE paUInt8
 typedef unsigned char SAMPLE;
 #define SAMPLE_SILENCE (128)
 #define PRINTF_S_FORMAT "%d"
 #endif

///MAKING SAMPLES READY
/*******************************************************************/
void dft(float inreal[],float outreal[],int n) {
int k,t;
    for ( k = 0; k < n; k++) {  // For each output element
        float sumreal = 0;
        for ( t = 0; t < n; t++) {  // For each input element
            double angle = 2 * PI * t * k / n;
            sumreal +=  inreal[t] * cos(angle);
        }
        outreal[k] = sumreal;
    }
}
void dodft(float* input,float* output,int len,int framesize){
	int i=0;
	while(i<len){
		dft(input+i,output+i,framesize);
		i+=882;
	}
	/*if(remain>0){
		dft(input,output,framesize);
		dodft(input+882,output+882,remain-882,882);
	}*/
}
/*******************************************************************/
int get_dft_input(float *samples,float* dftin,int num){
	FILE * pointer=fopen("voices.txt","w+");
	int k,index=0,i=0;
	for(k=0;k<num;++k){
		if((samples[k]>0.05) || (samples[k]<-0.05)){
			dftin[index]=samples[k];
			fprintf(pointer,"%f\n",samples[k]);
			++index;
			i=k;
		}
	}
	while(index%882!=0){
		dftin[index]=samples[i];
		fprintf(pointer,"%f\n",samples[i]);
		++index;
		++i;
	}
	fclose(pointer);
	return index;
}
/*****************************************************************/
void write_dft_result(float* dftout,int len,FILE* file){
	int j=0;
	file=fopen("dftresult.txt","w+");
	for(j=0;j<len;++j){
		fprintf(file,"%f\n",dftout[j]);
	}
	fclose(file);
}
/*****************************************************************/
#if USE_FRAME
void ham(float* dftin,int num){
	int i=0;
	int j=0;
	while(i<=num-882){
		j=0;
		while(j<882){
			dftin[i+j]=dftin[j+i]*(0.53836-(0.46164*cos((2* PI * j) /881)));
			++j;
		}
		i+=882;
	}
}
#else
void ham(float* dftin,int num){
	int i = 0;
	while(i<num){
		dftin[i]=dftin[i]*(0.53836-(0.46164*cos((2* PI * i) /(num-1))));
		++i;
	}
}
#endif
/*****************************************************************/
void emphasis(float* array,int len){
	float* backup=(float*)malloc(sizeof(float)*len);
	float a=0.95;
	backup[0]=array[0];
	int i=1;
	for(;i<len;++i){
		backup[i]=array[i]-a*array[i-1];
	}
	for(i=0;i<len;++i){
		array[i]=backup[i];
	}
	free(backup);
}
/*****************************************************************/
void copy(float*input,float* output,int istart,int ostart){
	int i=0;
	for(;i<882;++i){
		output[i+ostart]=input[i+istart];
	}
}
void frame(float* input,float*output,int len){
	int i=0;
	int j=0;
	while(i<=len-882){
		copy(input,output,i,j);
		i+=441;
		j+=882;
	}
}

/*****************************************************************/
int get_dft_input(float[],float[],int);
void ham(float[],int);
void write_dft_result(float[],int,FILE*);
void dft(float[],float[],int);
void emphasis(float[],int);
void frame(float[],float[],int);
void copy(float[],float[],int,int);

#if IS_UBUNTU
 /*
  * $Id$
  *
  * This program uses the PortAudio Portable Audio Library.
  * For more information see: http://www.portaudio.com
  * Copyright (c) 1999-2000 Ross Bencina and Phil Burk
  *
  * Permission is hereby granted, free of charge, to any person obtaining
  * a copy of this software and associated documentation files
  * (the "Software"), to deal in the Software without restriction,
  * including without limitation the rights to use, copy, modify, merge,
  * publish, distribute, sublicense, and/or sell copies of the Software,
  * and to permit persons to whom the Software is furnished to do so,
  * subject to the following conditions:
  *
  * The above copyright notice and this permission notice shall be
  * included in all copies or substantial portions of the Software.
  *
  * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
  * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
  * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
  * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
  * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
  */

 /*
34  * The text above constitutes the entire PortAudio license; however,
35  * the PortAudio community also makes the following non-binding requests:
36  *
37  * Any person wishing to distribute modifications to the Software is
38  * requested to send the modifications to the original developer so that
39  * they can be incorporated into the canonical version. It is also
40  * requested that these non-binding requests be included along with the
41  * license above.
42  */


 #include "portaudio.h"

 typedef struct
 {
  int frameIndex; /* Index into sample array. */
  int maxFrameIndex;
  SAMPLE *recordedSamples;
 }
 paTestData;

 /* This routine will be called by the PortAudio engine when audio is needed.
90 ** It may be called at interrupt level on some machines so don't do anything
91 ** that could mess up the system like calling malloc() or free().
92 */
 static int recordCallback( const void *inputBuffer, void *outputBuffer,
  unsigned long framesPerBuffer,
  const PaStreamCallbackTimeInfo* timeInfo,
  PaStreamCallbackFlags statusFlags,
  void *userData )
{
 paTestData *data = (paTestData*)userData;
  const SAMPLE *rptr = (const SAMPLE*)inputBuffer;
  SAMPLE *wptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
  long framesToCalc;
  long i;
  int finished;
  unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;

  (void) outputBuffer; /* Prevent unused variable warnings. */
  (void) timeInfo;
  (void) statusFlags;
  (void) userData;

  if( framesLeft < framesPerBuffer )
  {
  framesToCalc = framesLeft;
  finished = paComplete;
  }
  else
  {
  framesToCalc = framesPerBuffer;
  finished = paContinue;
  }

  if( inputBuffer == NULL )
  {
  for( i=0; i<framesToCalc; i++ )
  {
  *wptr++ = SAMPLE_SILENCE; /* left */
  if( NUM_CHANNELS == 2 ) *wptr++ = SAMPLE_SILENCE; /* right */
  }
  }
  else
  {
  for( i=0; i<framesToCalc; i++ )
 {
  *wptr++ = *rptr++; /* left */
  if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++; /* right */
  }
  }
  data->frameIndex += framesToCalc;
  return finished;
 }

 SAMPLE* record(int numSamples,int totalFrames){

  PaStreamParameters inputParameters,
  outputParameters;
  PaStream* stream;
  PaError err = paNoError;
  paTestData data;
  int i;
  int numBytes;
  SAMPLE max, val;
  double average;

  data.maxFrameIndex = totalFrames; /* Record for a few seconds. */
  data.frameIndex = 0;
  numBytes = numSamples * sizeof(SAMPLE);
  data.recordedSamples = (SAMPLE *) malloc( numBytes ); /* From now on, recordedSamples is initialised. */
  if( data.recordedSamples == NULL )
  {
  printf("Could not allocate record array.\n");
  goto done;
  }
  for( i=0; i<numSamples; i++ ) data.recordedSamples[i] = 0;

  err = Pa_Initialize();
  if( err != paNoError ) goto done;

  inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
  if (inputParameters.device == paNoDevice) {
  fprintf(stderr,"Error: No default input device.\n");
  goto done;
  }
  inputParameters.channelCount = 2; /* stereo input */
  inputParameters.sampleFormat = PA_SAMPLE_TYPE;
  inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
  inputParameters.hostApiSpecificStreamInfo = NULL;

  /* Record some audio. -------------------------------------------- */
  err = Pa_OpenStream(
  &stream,
  &inputParameters,
  NULL, /* &outputParameters, */
  SAMPLE_RATE,
  FRAMES_PER_BUFFER,
  paClipOff, /* we won't output out of range samples so don't bother clipping them */
  recordCallback,
  &data );
  if( err != paNoError ) goto done;

  err = Pa_StartStream( stream );
  if( err != paNoError ) goto done;
	#if IS_UBUNTU
	system("clear");
	#else
	system("cls");
	#endif
  printf("=== Now recording!! Please speak into the microphone. ===\n"); fflush(stdout);

  while( ( err = Pa_IsStreamActive( stream ) ) == 1 );
  if( err < 0 ) goto done;
  err = Pa_CloseStream( stream );
 done:
  Pa_Terminate();
  if( err != paNoError )
  {
  fprintf( stderr, "An error occured while using the portaudio stream\n" );
  fprintf( stderr, "Error number: %d\n", err );
  fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
  }
  return data.recordedSamples;
}

static int playCallback( const void *inputBuffer, void *outputBuffer,
                            unsigned long framesPerBuffer,
                            const PaStreamCallbackTimeInfo* timeInfo,
                            PaStreamCallbackFlags statusFlags,
                            void *userData )
   {
       paTestData *data = (paTestData*)userData;
       SAMPLE *rptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
       SAMPLE *wptr = (SAMPLE*)outputBuffer;
       unsigned int i;
       int finished;
       unsigned int framesLeft = data->maxFrameIndex - data->frameIndex;

       (void) inputBuffer; /* Prevent unused variable warnings. */
       (void) timeInfo;
       (void) statusFlags;
       (void) userData;

       if( framesLeft < framesPerBuffer )
       {
           /* final buffer... */
           for( i=0; i<framesLeft; i++ )
           {
               *wptr++ = *rptr++;  /* left */
               if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
           }
           for( ; i<framesPerBuffer; i++ )
           {
               *wptr++ = 0;  /* left */
               if( NUM_CHANNELS == 2 ) *wptr++ = 0;  /* right */
           }
           data->frameIndex += framesLeft;
           finished = paComplete;
       }
       else
       {
           for( i=0; i<framesPerBuffer; i++ )
           {
               *wptr++ = *rptr++;  /* left */
               if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
           }
           data->frameIndex += framesPerBuffer;
           finished = paContinue;
       }
       return finished;
   }

void playSample(SAMPLE *Sample){
  PaStreamParameters inputParameters,
  outputParameters;
  PaStream* stream;
  PaError err = paNoError;
  paTestData data;
  int i;
  int numBytes;
  SAMPLE max, val;
  double average;
    data.maxFrameIndex = TOTALFRAMES; /* Record for a few seconds. */
  data.frameIndex = 0;
  numBytes = NUMSAMPLES * sizeof(SAMPLE);
  data.recordedSamples = Sample;
    data.frameIndex = 0;

	err = Pa_Initialize();
  if( err != paNoError ) goto done;

    outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
    if (outputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default output device.\n");
        goto done;
    }
    outputParameters.channelCount = 2;                     /* stereo output */
    outputParameters.sampleFormat =  PA_SAMPLE_TYPE;
    outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;
	#if IS_UBUNTU
	system("clear");
	#else
	system("cls");
	#endif
    printf("\n=== Now playing back. ===\n"); fflush(stdout);
    err = Pa_OpenStream(
                &stream,
                NULL, /* no input */
                &outputParameters,
                SAMPLE_RATE,
                FRAMES_PER_BUFFER,
                paClipOff,      /* we won't output out of range samples so don't bother clipping them */
                playCallback,
                &data );
    if( err != paNoError ) goto done;
       if( stream )
       {
           err = Pa_StartStream( stream );
           if( err != paNoError ) goto done;

           printf("Waiting for playback to finish.\n"); fflush(stdout);

           while( ( err = Pa_IsStreamActive( stream ) ) == 1 ) Pa_Sleep(100);
           if( err < 0 ) goto done;

           err = Pa_CloseStream( stream );
           if( err != paNoError ) goto done;

           printf("Done.\n"); fflush(stdout);
       }
done:
  Pa_Terminate();
  if( err != paNoError )
  {
  fprintf( stderr, "An error occured while using the portaudio stream\n" );
  fprintf( stderr, "Error number: %d\n", err );
  fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
  }
}
#endif

int mod(int itot,int idiv){
    int it,iret;
    it=itot/idiv;
    iret = itot - idiv*it;
    return iret;
}

/*void dft(float datai[],float datao[],int n,int m1,int m2){
    double pi=3.141592653589793,arg,ct,st;
    int i,it,m,mt,imt,itex,zst;
    mt=0;
    for(m=m1; m<=m2; m+=1){
        datao[mt]=0;
        datao[mt+1]=0;
        it=0;
        for(i=-n/2; i<n/2; i+=1){
            imt=mod(i*m,n);
            arg=2*pi*imt/n;
            ct=cos(arg);
            st=sin(arg);
            datao[mt]+=datai[it]*ct-datai[it+1]*st;
            datao[mt+1]+=datai[it+1]*ct+datai[it]*st;
            it+=2;
        }
        mt+=2;
    }
    return;
}*/

int readSample(float SampleArr[], char filename[]){
    FILE *f;
	float s;
	f = fopen(filename,"r");
	int i = 0;
	if(f){
		while(fscanf(f, "%f", &s) != EOF)
			SampleArr[i++] = s;
		fclose(f);
		return 1;
	}
	return 0;
}


float dot_product(float v[], float u[], int n){
    int i;
    float result = 0;
    for (i = 0; i < n; i++)
    {
        result += v[i]*u[i];
    }
    return result;
}

float Vsize(float v[], int n){
    int i;
    float result = 0;
    for(i = 0; i < n; i++)
        result += pow(v[i],2);
    return sqrt(result);
}

int cutSample(float s1[], int size){
    int i, j = 0;
	for(i = 0; i < size; i++){
		if(fabs(s1[i]) > 0.15)
			s1[j++] = s1[i];	
	}
	return j;
}


int cutSamples(float s1[], float s2[], int size){
    int s1start, s2start, s1end, s2end;

    //sample1 starting point
    s1start = 5000;
    while(fabs(s1[s1start]) < 0.04)
        s1start++;
    while(fabs(s1[s1start]) > 0.002)
        s1start--;
    s1start++;

    //sample1 ending point
    s1end = size - 1;
    while(fabs(s1[s1end]) < 0.04)
        s1end--;
    while(fabs(s1[s1end]) > 0.002)
        s1end++;
    s1end--;

    //sample2 starting point
    s2start = 5000;
    while(fabs(s2[s2start]) < 0.04)
        s2start++;
    while(fabs(s2[s2start]) > 0.002)
        s2start--;
    s2start++;

    s2end = s2start + (s1end - s1start);

    int i;
    for(i = 0; i < s1end - s1start; i++){
        s1[i] = s1[i + s1start];
        s2[i] = s2[i + s2start];
    }
    #if TEST_PRINTS
    printf("s1start: %d\ns1end: %d\ns2start: %d\ns2end: %d\nlen %d\n\n", s1start, s1end, s2start, s2end, s1end - s1start);
    #endif
    return s1end - s1start;
}

float sumabs(float array[], int start, int len){
    int i;
    float sum = 0;
    for(i = start; i < len ;i++)
        sum = sum + fabs(array[i]) ;
    return sum ;
}

float sumabsdtft(float array[], int start, int len){
    int i;
    float sum = 0;
    for(i = start; i < len ;i += 2)
        sum = sum + fabs(array[i] + array[i + 1]);
    return sum ;
}

int writeSample(float SampleArr[], char filename[], int numSamples){
    FILE * f;
	float s;
	f = fopen(filename,"w+");
	int i;
	if(f){
		for(i = 0; i < numSamples; i++){
            fprintf(f, "%f\n", SampleArr[i]);
		}
		fclose(f);
		return 1;
	}
	return 0;
}

void GetSample(){
	#if IS_UBUNTU
	system("clear");
	#else
	system("cls");
	#endif

    //get sample name
    puts("Enter Sample Name (or leave empty to return)");
    char SampleName[100], c;
    int i = 0;
    while((c = getchar()) != '\n')
        SampleName[i++] = c;
    if(i == 0)
        return;
    SampleName[i] = '\0';
    while(strcmp(SampleName, "Samples") == 0 || strcmp(SampleName, "Samplestemp") == 0){
        puts("This is invalid. Enter a valid Sample Name (or leave empty to return)");
        i = 0;
        while((c = getchar()) != '\n')
            SampleName[i++] = c;
        if(i == 0)
            return;
        SampleName[i] = '.';
    }

    //add .txt
    SampleName[i++] = '.';
    SampleName[i++] = 't';
    SampleName[i++] = 'x';
    SampleName[i++] = 't';
    SampleName[i] = 0;

    FILE *f;

    //check if sample exists
    char SampleExists = 0;
    f = fopen("Samples.txt","r");
    if(f){

        //checking if sample exists
        char check[100];
        while(fscanf(f, "%s", check) != EOF){
            if(strcmp(check,SampleName) == 0){
                SampleExists = 1;
                printf("Sample Exists. Overwrite?(Y/N) ");
                c = getch();
                while(c != 'n' && c != 'N' && c != 'Y' && c != 'y')
                    c = getch();
                putchar('\n');
                if(c == 'n' || c == 'N')
                    return;
            }
        }
        fclose(f);
    }

    #if IS_UBUNTU
    int totalFrames = NUM_SECONDS * SAMPLE_RATE;
    int numSamples = totalFrames * NUM_CHANNELS;
    puts("press any key to start recording...");
    getch();
	SAMPLE *Sample = record(numSamples ,totalFrames);
	writeSample(Sample, SampleName, numSamples);
	SampleName[strlen(SampleName) - 3] = 'd';
	SampleName[strlen(SampleName) - 2] = 'a';
	SampleName[strlen(SampleName) - 1] = 't';

    double *spectrum=(double*)malloc(sizeof(double)*400000);
    FILE *sampleFile;
    unsigned int coeff;
    double mfcc_result;
	
	puts("processing...");
	#if USE_FRAME
	float *dftout=(float *)malloc(sizeof(float)*NUMSAMPLES*3);
	float *framed=(float *)malloc(sizeof(float)*NUMSAMPLES*3);
	int sizeAfterFraming=2*NUMSAMPLES - 882;
	emphasis(Sample,NUMSAMPLES);
	frame(Sample,framed,NUMSAMPLES);
	ham(framed,sizeAfterFraming);
	dodft(framed,dftout,sizeAfterFraming,882);
	free(framed);
    int j;
	for(j = 0; j < sizeAfterFraming; j++){
		spectrum[j] = dftout[j];
	}
	#else
	int cuttedSize = cutSample(Sample, NUMSAMPLES);
	emphasis(Sample,cuttedSize);
	ham(Sample,cuttedSize);
	float *dftout=(float *)malloc(sizeof(float)*NUMSAMPLES);
	dft(Sample,dftout,cuttedSize);
    int j;
	for(j = 0; j < cuttedSize; j++){
		spectrum[j] = dftout[j];
	}
	#endif
		double sum=0;
		puts("writing MFCC features");
		FILE* file1=fopen(SampleName,"w+");
		for(coeff = 1; coeff <= COEFFFS; coeff++)
		{
			printf("getting Coefficient %d/%d\n", coeff, COEFFFS);
			#if USE_DTW
			mfcc_result = GetCoefficient(spectrum, 44100, 40, sizeAfterFraming, coeff);
			#else
			mfcc_result = GetCoefficient(spectrum, 44100, 1000, 512, coeff);
			#endif
			fprintf(file1,"%f\n",fabs(mfcc_result));
		}
		fclose(file1);
		free(Sample);
		free(spectrum);
		free(dftout);
		#endif
		f = fopen("Samples.txt","a+");
        if(f){
            if(SampleExists == 0){
				SampleName[strlen(SampleName) - 3] = 't';
				SampleName[strlen(SampleName) - 2] = 'x';
				SampleName[strlen(SampleName) - 1] = 't';
                fprintf(f,"%s\n",SampleName);
			}
            fclose(f);
            puts("Sample saved.");
        }
	puts("Press any key to continue...");
    getch();
}

int IsEqual(SAMPLE Sample1In[], SAMPLE Sample2In[]){ //Sum of Absolutes of DFT Coefficients
    /*SAMPLE Sample1Cut[NUMSAMPLES], Sample1Out[150];
    SAMPLE Sample2Cut[NUMSAMPLES], Sample2Out[150];

    int cuttedSize = cutSamples(Sample1In, Sample2In, Sample1Cut, Sample2Cut, NUMSAMPLES);
    int i,j;
    float num1, num2;


    float sum1 = sumabs(Sample1Out, 0, 122);
    float sum2 = sumabs(Sample2Out, 0, 122);
    float sum3 = sumabs(Sample1Cut, 0, cuttedSize);
    float sum4 = sumabs(Sample2Cut, 0, cuttedSize);

    float sumsample1 = 0 , sumsample2 = 0;


    dft(Sample1Cut, Sample1Out, cuttedSize / 2, -30, 30);
    dft(Sample2Cut, Sample2Out, cuttedSize / 2, -30, 30);
    #if TEST_PRINTS
    printf("before: %f %f %f\n", sum4, sum3, fabs(sum4 - sum3));
    printf("dft %f %f %f\n", sumabs(Sample1Out, 0, cuttedSize), sumabs(Sample2Out, 0, cuttedSize), fabs(sumabs(Sample1Out, 0, cuttedSize) - sumabs(Sample2Out, 0, cuttedSize)));
    putchar('\n');
    #endif

        sumsample1 = sumsample2 = 0;
        for(j = 0 ; j < cuttedSize ; j++){
            sumsample1 += (Sample1Cut[j] * Sample1Cut[j]);
            sumsample2 += (Sample2Cut[j] * Sample2Cut[j]);
        }
        for(j = 0 ; j < cuttedSize ; j++){
            Sample1Cut[j] = Sample1Cut[j] / sqrt(sumsample1);
            Sample2Cut[j] = Sample2Cut[j] / sqrt(sumsample2);
        }

        num1 = sumabs(Sample1Cut, 0, cuttedSize);
        num2 = sumabs(Sample2Cut, 0, cuttedSize);

    dft(Sample1Cut, Sample1Out, cuttedSize / 2, -30, 30);
    dft(Sample2Cut, Sample2Out, cuttedSize / 2, -30, 30);
    #if TEST_PRINTS
    printf("levelamirerfan: %f %f %f\n", num1, num2, fabs(num2 - num1));
    printf("dft %f %f %f\n\n", sumabs(Sample1Out, 0, 122), sumabs(Sample2Out, 0, 122), fabs(sumabs(Sample1Out, 0, 122) - sumabs(Sample2Out, 0, 122)));
    #endif
    if(fabs(sumabs(Sample1Out, 0, 122) - sumabs(Sample2Out, 0, 122)) < 10)*/
	#if TEST_PRINTS
    printf("%f", dot_product(Sample1In, Sample2In, COEFFFS) / (Vsize(Sample1In, COEFFFS) * Vsize(Sample2In, COEFFFS)));
	#endif
    if(dot_product(Sample1In, Sample2In, COEFFFS) / (Vsize(Sample1In, COEFFFS) * Vsize(Sample2In, COEFFFS)) > 0.1)
        return 1;
    return 0;
}

int SampleList(){
    FILE *f;
    f = fopen("Samples.txt","r");
    int samplesCount = 0;
    if(f){

		#if IS_UBUNTU
		system("clear");
		#else
		system("cls");
		#endif

        //writing samples on screen
        char temp[100];
        while(fscanf(f, "%s", temp) != EOF){
            temp[strlen(temp) - 4] = '\0'; //removing .txt
            samplesCount++;
            printf("%d. %s\n", samplesCount, temp);
        }
        fclose(f);
    }
    return samplesCount;
}

void ShowSamples(){
    int samplesCount = SampleList();
    if(samplesCount == 0)
        puts("No samples found to show.");
	else{
		int Sample1Choice = 0;
        while(Sample1Choice < 1 || Sample1Choice > samplesCount){
            printf("enter a number for sample to play: (or leave empty to return)");
            Sample1Choice = 0;
            char c = getchar();
            if(c == '\n')
                return;
            while(c >= '0' && c <= '9'){
                Sample1Choice = Sample1Choice * 10 + (c - '0');
                c = getchar();
            }
        }

		FILE *f;
        f = fopen("Samples.txt","r");
        //read samples
        char Sample1Filename[100];
        fseek(f,0,SEEK_SET);
        int i = 0;
        while(i < Sample1Choice){
            i++;
            fscanf(f, "%s", Sample1Filename);
        }

        SAMPLE Sample1[NUMSAMPLES];
        if(readSample(Sample1, Sample1Filename)){
			#if IS_UBUNTU
			playSample(Sample1);
			#endif // IS_UBUNTU
        }
        else{
            puts("Could not read samples");
        }
        fclose(f);
	}
	puts("Press any key to continue...");
    getch();
}

void CompareSamplesList(){
    int samplesCount = SampleList();
    if(samplesCount < 2)
        puts("Not enough samples found to compare.");
    else{

        //let user choose two samples
		char c;
        int Sample1Choice = 0, Sample2Choice = 0;
        while(Sample1Choice < 1 || Sample1Choice > samplesCount){
            printf("enter a number for sample1: ");
			Sample1Choice = 0;
			while((c = getchar()) >= '0' && c <= '9')
				Sample1Choice = Sample1Choice * 10 + c - '0';
        }
        while(Sample2Choice < 1 || Sample2Choice > samplesCount){
            printf("enter a number for sample2: ");
			Sample2Choice = 0;
			while((c = getchar()) >= '0' && c <= '9')
				Sample2Choice = Sample2Choice * 10 + c - '0';
        }
        FILE *f;
        f = fopen("Samples.txt","r");
        //read samples
        char Sample1Filename[100], Sample2Filename[100];
        fseek(f,0,SEEK_SET);
        int i = 0;
        while(i < Sample1Choice){
            i++;
            fscanf(f, "%s", Sample1Filename);
        }
        fseek(f,0,SEEK_SET);
        i = 0;
        while(i < Sample2Choice){
            i++;
            fscanf(f, "%s", Sample2Filename);
        }
        fclose(f);

		Sample1Filename[strlen(Sample1Filename) - 3] = 'd';
		Sample1Filename[strlen(Sample1Filename) - 2] = 'a';
		Sample1Filename[strlen(Sample1Filename) - 1] = 't';

		Sample2Filename[strlen(Sample2Filename) - 3] = 'd';
		Sample2Filename[strlen(Sample2Filename) - 2] = 'a';
		Sample2Filename[strlen(Sample2Filename) - 1] = 't';
		
		#if USE_DTW
		float dtwanswer = dtw(Sample1Filename, Sample2Filename, "output.txt", 12, 12, 1);
		#if TEST_PRINTS
		printf("%f\n", dtwanswer);
		#endif
		if(dtwanswer < DTW_PASS_VALUE)
			puts("PASS");
		else
			puts("FAIL");
		#else
        float Sample1[COEFFFS], Sample2[COEFFFS];
        if(readSample(Sample1, Sample1Filename) && readSample(Sample2, Sample2Filename)){
            puts("Trying Method1...");
            if(IsEqual(Sample1, Sample2))
                puts("PASS");
            else
                puts("FAIL");
        }
        else{
            puts("Could not read samples");
        }
		#endif
    }
	puts("Press any key to continue...");
    getch();
}

void DeleteSample(){
    int samplesCount = SampleList();
    if(samplesCount == 0)
        puts("No Sample found to delete.");
    else{

        //let user choose a sample
        int Sample1Choice = 0;
        while(Sample1Choice < 1 || Sample1Choice > samplesCount){
            printf("enter a number for sample: (or leave empty to return)");
            Sample1Choice = 0;
            char c = getchar();
            if(c == '\n')
                return;
            while(c >= '0' && c <= '9'){
                Sample1Choice = Sample1Choice * 10 + (c - '0');
                c = getchar();
            }
        }

        FILE *f1, *f2;
        f1 = fopen("Samples.txt","r");
        f2 = fopen("Samplestemp.txt","w");

        //make new list
        char Sample1Filename[100];
        fseek(f1,0,SEEK_SET);
        int i = 0;
        while(i < samplesCount){
            fscanf(f1, "%s", Sample1Filename);
            i++;
            if(i != Sample1Choice)
                fprintf(f2, "%s\n", Sample1Filename);
            else{
                remove(Sample1Filename);
				Sample1Filename[strlen(Sample1Filename) - 3] = 'd';
				Sample1Filename[strlen(Sample1Filename) - 2] = 'a';
				Sample1Filename[strlen(Sample1Filename) - 1] = 't';
                remove(Sample1Filename);
			}
        }
        fclose(f1);
        fclose(f2);
        remove("Samples.txt");
        rename("Samplestemp.txt", "Samples.txt");
    }
	puts("Done.\nPress any key to continue...");
    getch();
}

void FindMatches(){
	char c;
    int samplesCount = SampleList();
    if(samplesCount < 2)
        puts("Not enough samples found to compare.");
    else{

        //let user choose two samples
        int Sample1Choice = 0, Sample2Choice = 0;
        while(Sample1Choice < 1 || Sample1Choice > samplesCount){
            printf("enter a number for sample1: ");
			Sample1Choice = 0;
			while((c = getchar()) >= '0' && c <= '9')
				Sample1Choice = Sample1Choice * 10 + c - '0';
        }

        FILE *f;
        f = fopen("Samples.txt","r");
        //read samples
        char Sample1Filename[100], Sample2Filename[100];
        fseek(f,0,SEEK_SET);
        int i = 0;
        while(i < Sample1Choice){
            i++;
            fscanf(f, "%s", Sample1Filename);
        }
        fseek(f,0,SEEK_SET);
        i = 0;
        SAMPLE Sample1[COEFFFS], Sample2[COEFFFS];
		Sample1Filename[strlen(Sample1Filename) - 3] = 'd';
		Sample1Filename[strlen(Sample1Filename) - 2] = 'a';
		Sample1Filename[strlen(Sample1Filename) - 1] = 't';
        readSample(Sample1, Sample1Filename);
        while(i < samplesCount){
            i++;
            fscanf(f, "%s", Sample2Filename);
			Sample2Filename[strlen(Sample2Filename) - 3] = 'd';
			Sample2Filename[strlen(Sample2Filename) - 2] = 'a';
			Sample2Filename[strlen(Sample2Filename) - 1] = 't';
			#if USE_DTW
			float dtwanswer = dtw(Sample1Filename, Sample2Filename, "output.txt", 12, 12, 1);
			#if TEST_PRINTS
			printf("%f\n", dtwanswer);
			#endif
			if(dtwanswer < DTW_PASS_VALUE){
				Sample2Filename[strlen(Sample2Filename) - 4] = '\0';
                puts(Sample2Filename);			
			}
			#else
            if(readSample(Sample2, Sample2Filename)){
                if(IsEqual(Sample1, Sample2)){
					Sample2Filename[strlen(Sample2Filename) - 4] = '\0';
                    puts(Sample2Filename);
				}
            }
            else{
                printf("Could not read sample %s", Sample2Filename);
            }
			#endif
        }
        fclose(f);
    }
	puts("Press any key to continue...");
    getch();
}

int main(){
    int choice = 0;
    while(choice != '6'){
		#if IS_UBUNTU
		system("clear");
		#else
		system("cls");
		#endif
        puts("Welcome\nChoose a number:\n1.Get sample\n2.Show/Play Samples\n3.Compare samples\n4.delete sample\n5.find matches\n6.exit");
        choice = getch();
        switch(choice){
            case '1':
                GetSample();
                break;
            case '2':
                ShowSamples();
                break;
            case '3':
                CompareSamplesList();
                break;
            case '4':
                DeleteSample();
                break;
            case '5':
                FindMatches();
                break;
            }
    }
    return 0;
}
