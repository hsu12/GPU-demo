#include <opencv2/core/core.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/opencv.hpp>
//#include <opencv2/gpu/gpu.hpp>
//#include <opencv2/cudaarithm.hpp>
#include "cudaarithm.hpp"
#include <stdio.h>
#include <iostream>
#include <ctime>
#include <sys/time.h>
using namespace std;
extern void median_filter_wrapper(const cv::Mat& input, cv::Mat& output);


int main()
{
	// Read input file (image)
	std::string imagePath = "/home/nvidia/CUDA/cuda-convolution/data/imagebilateral.png";
	cv::Mat input = cv::imread(imagePath,0);
	if(input.empty()) {
		std::cout<<"Could not load image. Check location and try again."<<std::endl;
		std::cin.get();
		return -1;
	}

	double running_sum = 0.0;
	int attempts = 10;

	cv::Size resize_size;
	resize_size.width = 960;
	resize_size.height = 800;
	cv::resize(input,input,resize_size);
	cv::Mat output_gpu(input.rows,input.cols,CV_8UC1);
	cv::Mat output_cpu(input.rows,input.cols,CV_8UC1);	

	// --------------- MEDIAN FILTER ---------------
        for (int ctr = 0; ctr < attempts; ctr++) {
		clock_t gpu_s = clock();
		median_filter_wrapper(input,output_gpu);
		cout<<ctr<<endl;
		clock_t gpu_e = clock();
		if (ctr > 0)
		    running_sum = running_sum + (double(gpu_e-gpu_s)*1000)/CLOCKS_PER_SEC;
	}
	std::cout << "GPU Accelerated Median Filter took " << running_sum/(attempts-1) << " ms.\n";	
	running_sum = 0.0;
	for (int ctr = 0; ctr < attempts; ctr++) {
		clock_t cpu_s = clock();
		cv::medianBlur(input,output_cpu,9);
		clock_t cpu_e = clock();
		if (ctr > 0)
		    running_sum = running_sum + (double(cpu_e-cpu_s)*1000)/CLOCKS_PER_SEC;
	}
	std::cout << "CPU Accelerated Median Filter took " << running_sum/(attempts-1) << " ms.\n";	

	cv::imshow("(MF) Output Image - GPU",output_gpu);
	cv::imwrite("gpu_median_result.png",output_gpu);
	cv::imshow("(MF) Output Image - CPU",output_cpu);
	cv::imwrite("cpu_median_result.png",output_cpu);
	cv::waitKey();

	

	return 0;
}
