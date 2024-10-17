#include "opencv2/opencv.hpp"
#include "resizeim/resizeim.h"
#include <iostream>

int main() {
  cv::Mat image = cv::imread("image.png");

  if (image.empty()) {
    std::cout << "Error loading image" << std::endl;
    return 1;
  }

  coder::array<unsigned char, 3U> im;
  coder::array<unsigned char, 3U> out;

  im.set_size(image.rows, image.cols, 3);

  // Convert to uint8_t array
  for (int i = 0; i < image.rows; i++) {
    for (int j = 0; j < image.cols; j++) {
      for (int k = 0; k < 3; k++) {
        im[i + image.rows * j + image.rows * image.cols * k] =
            image.at<cv::Vec3b>(i, j)[k];
      }
    }
  }

  resizeim(im, out);

  // Convert im to cv::Mat and display
  cv::Mat outImage(image.rows / 2, image.cols / 2, CV_8UC3);
  for (int i = 0; i < image.rows / 2; i++) {
    for (int j = 0; j < image.cols / 2; j++) {
      for (int k = 0; k < 3; k++) {
        outImage.at<cv::Vec3b>(i, j)[k] =
            out[i + (image.rows / 2) * j +
                (image.rows / 2) * (image.cols / 2) * k];
      }
    }
  }

  // Print the original image size
  std::cout << "Original Image Size: " << image.rows << "x" << image.cols
            << std::endl;

  // Print the resized image size
  std::cout << "Resized Image Size: " << outImage.rows << "x" << outImage.cols
            << std::endl;

  cv::imshow("Original Image", image);
  cv::imshow("Resized Image", outImage);
  cv::waitKey(0);

  return 0;
}
