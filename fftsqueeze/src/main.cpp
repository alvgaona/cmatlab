#include "fftsqueeze/fftsqueeze.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>

void generate_sin_wave(coder::array<double, 2U> &out) {
  const int num_samples = 1000;
  const double frequency = 10.0;
  const double sample_rate = 1000.0;
  const double amplitude = 1.0;

  out.set_size(1, num_samples);

  for (auto i = 0; i < num_samples; i++) {
    double t = static_cast<double>(i) / sample_rate;
    out[i] = amplitude * std::sin(2.0 * M_PI * frequency * t);
  }
}

int main() {
  coder::array<creal_T, 2U> s;
  coder::array<double, 2U> t;
  coder::array<double, 2U> x;
  coder::array<double, 1U> f;
  coder::array<double, 1U> window;

  window.set_size(128);

  std::vector<double> kaiser_window = {
      0.9403, 0.9421, 0.9440, 0.9457, 0.9475, 0.9492, 0.9509, 0.9526, 0.9542,
      0.9558, 0.9574, 0.9590, 0.9605, 0.9620, 0.9635, 0.9650, 0.9664, 0.9678,
      0.9691, 0.9705, 0.9718, 0.9730, 0.9743, 0.9755, 0.9767, 0.9778, 0.9790,
      0.9801, 0.9811, 0.9822, 0.9832, 0.9842, 0.9851, 0.9861, 0.9870, 0.9878,
      0.9887, 0.9895, 0.9902, 0.9910, 0.9917, 0.9924, 0.9931, 0.9937, 0.9943,
      0.9949, 0.9954, 0.9959, 0.9964, 0.9968, 0.9973, 0.9977, 0.9980, 0.9983,
      0.9986, 0.9989, 0.9992, 0.9994, 0.9995, 0.9997, 0.9998, 0.9999, 1.0000,
      1.0000, 1.0000, 1.0000, 0.9999, 0.9998, 0.9997, 0.9995, 0.9994, 0.9992,
      0.9989, 0.9986, 0.9983, 0.9980, 0.9977, 0.9973, 0.9968, 0.9964, 0.9959,
      0.9954, 0.9949, 0.9943, 0.9937, 0.9931, 0.9924, 0.9917, 0.9910, 0.9902,
      0.9895, 0.9887, 0.9878, 0.9870, 0.9861, 0.9851, 0.9842, 0.9832, 0.9822,
      0.9811, 0.9801, 0.9790, 0.9778, 0.9767, 0.9755, 0.9743, 0.9730, 0.9718,
      0.9705, 0.9691, 0.9678, 0.9664, 0.9650, 0.9635, 0.9620, 0.9605, 0.9590,
      0.9574, 0.9558, 0.9542, 0.9526, 0.9509, 0.9492, 0.9475, 0.9457, 0.9440,
      0.9421, 0.9403,
  };

  for (int i = 0; i < 128; ++i) {
    window[i] = kaiser_window[i];
  }

  generate_sin_wave(x);

  // Call the entry-point 'fftsqueeze'.
  fftsqueeze(x, 1000, window, s, f, t);
  // Convert the complex spectrogram to magnitude
  cv::Mat spectrogram(s.size(0), s.size(1), CV_64F);
  for (int i = 0; i < s.size(0); ++i) {
    for (int j = 0; j < s.size(1); ++j) {
      const creal_T value = s[i + j * s.size(0)];
      const std::complex<double> z(value.re, value.im);
      spectrogram.at<double>(i, j) = std::abs(z);
    }
  }

  // Normalize the spectrogram for better visualization
  cv::normalize(spectrogram, spectrogram, 0, 255, cv::NORM_MINMAX, CV_8U);

  // Apply a color map
  cv::Mat colorSpectrogram;
  cv::applyColorMap(spectrogram, colorSpectrogram, cv::COLORMAP_JET);

  // Flip the image vertically because OpenCV's origin is at the top-left
  cv::flip(colorSpectrogram, colorSpectrogram, 0);

  // Resize the spectrogram to 800x600
  cv::resize(colorSpectrogram, colorSpectrogram, cv::Size(800, 600), 0, 0,
             cv::INTER_LINEAR);

  // Display the spectrogram
  cv::imshow("Spectrogram", colorSpectrogram);
  cv::waitKey(0);

  return 0;
}
