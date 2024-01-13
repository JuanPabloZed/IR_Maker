        /// <summary>
        /// Transforms input filter coefficients to minimum phase coefficients.
        /// </summary>
        /// <param name="input">Input filter coefficients.</param>
        /// <returns>Minimum phase transformed filter coefficients.</returns>
        public static double[] Transform(double[] input)
        {
            int n = input.Length;
            int power = (int)Math.Ceiling(Math.Log(n, 2d)) + 3; // Larger fft size instead double min required to prevent aliasing...
            int fftSize = (int)Math.Pow(2, power);
            int halfFttSize = fftSize / 2;
            int halfFftSizeP1 = halfFttSize + 1;
            double[] kernel = new double[fftSize];
            double[] kernelCss = new double[fftSize + 2];
            double[] magnitudes = new double[halfFftSizeP1];
            Array.Copy(input, kernel, n);
            // Double precision FFT transformer.
            FFTDouble fft = new FFTDouble(fftSize);

            // Forward transformation. Result is halfFttSize+1 real-imaginary pairs.
            fft.FftForward(kernel, 0, kernelCss, 0);

            int offset=0, offsetP1;
            double real, imaginary, magnitude;

            for (int i = 0; i < halfFftSizeP1; i++)
            {
                offsetP1 = offset + 1;
                real = kernelCss[offset];
                imaginary = kernelCss[offsetP1];
                magnitude = Math.Sqrt(real * real + imaginary * imaginary);
                magnitudes[i] = magnitude;
                kernelCss[offset] = Math.Log(magnitude);
                kernelCss[offsetP1] = 0d;
                offset += 2;
            }

            // Backward transformation.
            fft.FftBackward(kernelCss, 0, kernel, 0);

            for (int i = halfFttSize; i < fftSize; i++)
                kernel[i] *= -1d;

            fft.FftForward(kernel, 0, kernelCss, 0);

            offset = 0;

            for (int i = 0; i < halfFftSizeP1; i++)
            {
                offsetP1 = offset + 1;
                imaginary = kernelCss[offsetP1];
                real = Math.Cos(imaginary) * magnitudes[i];
                imaginary = Math.Sin(imaginary) * magnitudes[i];
                kernelCss[offset] = real;
                kernelCss[offsetP1] = imaginary;
                offset += 2;
            }

            fft.FftBackward(kernelCss, 0, kernel, 0);
            fft.Free();
            fft = null;
            double[] result = new double[n];
            Array.Copy(kernel, result, n);
            return result;
        }