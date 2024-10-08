using System;
using System.Diagnostics;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;

namespace DigitalMusicAnalysis
{
    public class timefreq
    {
        public float[][] timeFreqData;
        public int wSamp;
        public Complex[] twiddles;
        public Stopwatch stopwatch3 = new Stopwatch();

        public timefreq(float[] x, int windowSamp)
        {
            int ii;
            double pi = 3.14159265;
            Complex i = Complex.ImaginaryOne;
            this.wSamp = windowSamp;
            twiddles = new Complex[wSamp];

            for (ii = 0; ii < wSamp; ii++)
            {
                double a = 2 * pi * ii / (double)wSamp;
                twiddles[ii] = Complex.Pow(Complex.Exp(-i), (float)a);
            }

            timeFreqData = new float[wSamp/2][];

            int nearest = (int)Math.Ceiling((double)x.Length / (double)wSamp);
            nearest = nearest * wSamp;

            Complex[] compX = new Complex[nearest];
            for (int kk = 0; kk < nearest; kk++)
            {
                if (kk < x.Length)
                {
                    compX[kk] = x[kk];
                }
                else
                {
                    compX[kk] = Complex.Zero;
                }
            }


            int cols = 2 * nearest /wSamp;

            for (int jj = 0; jj < wSamp / 2; jj++)
            {
                timeFreqData[jj] = new float[cols];
            }

            timeFreqData = stft(compX, wSamp); //11%
	
        }

        float[][] stft(Complex[] x, int wSamp)
        {
            int N = x.Length;
            float fftMax = 0;
            float[][] Y = new float[wSamp / 2][];

            stopwatch3.Start();

            // Pre-initialize the 2D array
            for (int ll = 0; ll < wSamp / 2; ll++)
            {
                Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            }

            Complex[] temp = new Complex[wSamp];
            Complex[] tempFFT = new Complex[wSamp];

            Parallel.For(0, (int)(2 * Math.Floor((double)N / (double)wSamp) - 1), new ParallelOptions { MaxDegreeOfParallelism = MainWindow.numThreads }, ii =>
            {
                Complex[] localTemp = new Complex[wSamp]; // Thread-local temp array

                // Populate temp in parallel
                for (int jj = 0; jj < wSamp; jj++)
                {
                    localTemp[jj] = x[ii * (wSamp / 2) + jj];
                }

                // Perform FFT sequentially (due to FFT recursion)
                Complex[] localTempFFT = fft(localTemp);

                float localFftMax = 0; // Thread-local fftMax to avoid race conditions

                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] = (float)Complex.Abs(localTempFFT[kk]);
                    if (Y[kk][ii] > localFftMax)
                    {
                        localFftMax = Y[kk][ii];
                    }
                }

                // Atomic max update using Interlocked.CompareExchange
                float initialMax, computedMax;
                do
                {
                    initialMax = fftMax;
                    computedMax = Math.Max(initialMax, localFftMax);
                } while (initialMax != Interlocked.CompareExchange(ref fftMax, computedMax, initialMax));
            });

            // Normalize Y array by fftMax
            Parallel.For(0, 2 * (int)Math.Floor((double)N / (double)wSamp) - 1, new ParallelOptions { MaxDegreeOfParallelism = MainWindow.numThreads }, ii =>
            {
                for (int kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] /= fftMax;
                }
            });

            stopwatch3.Stop();

            return Y;
        }


        Complex[] fft(Complex[] x)
        {
            int ii = 0;
            int kk = 0;
            int N = x.Length;

            Complex[] Y = new Complex[N];


            if (N == 1)
            {
                Y[0] = x[0];
            }
            else{

                Complex[] E = new Complex[N/2];
                Complex[] O = new Complex[N/2];
                Complex[] even = new Complex[N/2];
                Complex[] odd = new Complex[N/2];

                for (ii = 0; ii < N; ii++)
                {

                    if (ii % 2 == 0)
                    {
                        even[ii / 2] = x[ii];
                    }
                    if (ii % 2 == 1)
                    {
                        odd[(ii - 1) / 2] = x[ii];
                    }
                }

                E = fft(even);
                O = fft(odd);

                for (kk = 0; kk < N; kk++)
                {
                    Y[kk] = E[(kk % (N / 2))] + O[(kk % (N / 2))] * twiddles[kk * wSamp / N];
                }
            }

           return Y;
        }
        
    }
}
