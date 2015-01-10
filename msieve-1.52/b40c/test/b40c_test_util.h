/******************************************************************************
 * Copyright 2010 Duane Merrill
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. 
 * 
 * For more information, see our Google Code project site: 
 * http://code.google.com/p/back40computing/
 * 
 * Thanks!
 ******************************************************************************/

#pragma once

#if defined(_WIN32) || defined(_WIN64)
	#include <windows.h>
	#undef small			// Windows is terrible for polluting macro namespace
#else
	#include <sys/resource.h>
#endif



#include <stdio.h>
#include <math.h>
#include <float.h>

#include <map>
#include <string>
#include <sstream>
#include <iostream>

#include <b40c/util/random_bits.cuh>
#include <b40c/util/basic_utils.cuh>

namespace b40c {


/******************************************************************************
 * Command-line parsing
 ******************************************************************************/

class CommandLineArgs
{
protected:

	std::map<std::string, std::string> pairs;

public:

	// Constructor
	CommandLineArgs(int argc, char **argv)
	{
		using namespace std;

	    for (int i = 1; i < argc; i++)
	    {
	        string arg = argv[i];

	        if ((arg[0] != '-') || (arg[1] != '-')) {
	        	continue;
	        }

        	string::size_type pos;
		    string key, val;
	        if ((pos = arg.find( '=')) == string::npos) {
	        	key = string(arg, 2, arg.length() - 2);
	        	val = "";
	        } else {
	        	key = string(arg, 2, pos - 2);
	        	val = string(arg, pos + 1, arg.length() - 1);
	        }
        	pairs[key] = val;
	    }
	}

	bool CheckCmdLineFlag(const char* arg_name)
	{
		using namespace std;
		map<string, string>::iterator itr;
		if ((itr = pairs.find(arg_name)) != pairs.end()) {
			return true;
	    }
		return false;
	}

	template <typename T>
	void GetCmdLineArgument(const char *arg_name, T &val);

	int ParsedArgc()
	{
		return pairs.size();
	}
};

template <typename T>
void CommandLineArgs::GetCmdLineArgument(const char *arg_name, T &val)
{
	using namespace std;
	map<string, string>::iterator itr;
	if ((itr = pairs.find(arg_name)) != pairs.end()) {
		istringstream strstream(itr->second);
		strstream >> val;
    }
}

template <>
void CommandLineArgs::GetCmdLineArgument<char*>(const char* arg_name, char* &val)
{
	using namespace std;
	map<string, string>::iterator itr;
	if ((itr = pairs.find(arg_name)) != pairs.end()) {

		string s = itr->second;
		val = (char*) malloc(sizeof(char) * (s.length() + 1));
		strcpy(val, s.c_str());

	} else {
    	val = NULL;
	}
}





/******************************************************************************
 * Device initialization
 ******************************************************************************/

void DeviceInit(CommandLineArgs &args)
{
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	if (deviceCount == 0) {
		fprintf(stderr, "No devices supporting CUDA.\n");
		exit(1);
	}
	int dev = 0;
	args.GetCmdLineArgument("device", dev);
	if (dev < 0) {
		dev = 0;
	}
	if (dev > deviceCount - 1) {
		dev = deviceCount - 1;
	}
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, dev);
	if (deviceProp.major < 1) {
		fprintf(stderr, "Device does not support CUDA.\n");
		exit(1);
	}
    if (!args.CheckCmdLineFlag("quiet")) {
        printf("Using device %d: %s\n", dev, deviceProp.name);
    }

	cudaSetDevice(dev);
}




/******************************************************************************
 * Templated routines for printing keys/values to the console 
 ******************************************************************************/

template<typename T> 
void PrintValue(T val) {
	val.Print();
}

template<>
void PrintValue<char>(char val) {
	printf("%d", val);
}

template<>
void PrintValue<short>(short val) {
	printf("%d", val);
}

template<>
void PrintValue<int>(int val) {
	printf("%d", val);
}

template<>
void PrintValue<long>(long val) {
	printf("%ld", val);
}

template<>
void PrintValue<long long>(long long val) {
	printf("%lld", val);
}

template<>
void PrintValue<float>(float val) {
	printf("%f", val);
}

template<>
void PrintValue<double>(double val) {
	printf("%f", val);
}

template<>
void PrintValue<unsigned char>(unsigned char val) {
	printf("%u", val);
}

template<>
void PrintValue<unsigned short>(unsigned short val) {
	printf("%u", val);
}

template<>
void PrintValue<unsigned int>(unsigned int val) {
	printf("%u", val);
}

template<>
void PrintValue<unsigned long>(unsigned long val) {
	printf("%lu", val);
}

template<>
void PrintValue<unsigned long long>(unsigned long long val) {
	printf("%llu", val);
}



/******************************************************************************
 * Helper routines for list construction and validation 
 ******************************************************************************/

/**
 * Compares the equivalence of two arrays
 */
template <typename T, typename SizeT>
int CompareResults(T* computed, T* reference, SizeT len, bool verbose = true)
{
	for (SizeT i = 0; i < len; i++) {

		if (computed[i] != reference[i]) {
			printf("INCORRECT: [%lu]: ", (unsigned long) i);
			PrintValue<T>(computed[i]);
			printf(" != ");
			PrintValue<T>(reference[i]);

			if (verbose) {
				printf("\nresult[...");
				for (size_t j = (i >= 5) ? i - 5 : 0; (j < i + 5) && (j < len); j++) {
					PrintValue<T>(computed[j]);
					printf(", ");
				}
				printf("...]");
				printf("\nreference[...");
				for (size_t j = (i >= 5) ? i - 5 : 0; (j < i + 5) && (j < len); j++) {
					PrintValue<T>(reference[j]);
					printf(", ");
				}
				printf("...]");
			}

			return 1;
		}
	}

	printf("CORRECT");
	return 0;
}


/**
 * Verify the contents of a device array match those
 * of a host array
 */
template <typename T>
int CompareDeviceResults(
	T *h_reference,
	T *d_data,
	size_t num_elements,
	bool verbose = true,
	bool display_data = false)
{
	// Allocate array on host
	T *h_data = (T*) malloc(num_elements * sizeof(T));

	// Reduction data back
	cudaMemcpy(h_data, d_data, sizeof(T) * num_elements, cudaMemcpyDeviceToHost);

	// Display data
	if (display_data) {
		printf("Reference:\n");
		for (int i = 0; i < num_elements; i++) {
			PrintValue(h_reference[i]);
			printf(", ");
		}
		printf("\n\nData:\n");
		for (int i = 0; i < num_elements; i++) {
			PrintValue(h_data[i]);
			printf(", ");
		}
		printf("\n\n");
	}

	// Check
	int retval = CompareResults(h_data, h_reference, num_elements, verbose);

	// Cleanup
	if (h_data) free(h_data);

	return retval;
}

int CompareDeviceResults(
	util::NullType *h_reference,
	util::NullType *d_data,
	size_t num_elements,
	bool verbose = true,
	bool display_data = false)
{
	return 0;
}

/**
 * Verify the contents of a device array match those
 * of a host array
 */
template <typename T>
int CompareDeviceDeviceResults(
	T *d_reference,
	T *d_data,
	size_t num_elements,
	bool verbose = true,
	bool display_data = false)
{
	// Allocate array on host
	T *h_reference = (T*) malloc(num_elements * sizeof(T));
	T *h_data = (T*) malloc(num_elements * sizeof(T));

	// Reduction data back
	cudaMemcpy(h_reference, d_reference, sizeof(T) * num_elements, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_data, d_data, sizeof(T) * num_elements, cudaMemcpyDeviceToHost);

	// Display data
	if (display_data) {
		printf("Reference:\n");
		for (int i = 0; i < num_elements; i++) {
			PrintValue(h_reference[i]);
			printf(", ");
		}
		printf("\n\nData:\n");
		for (int i = 0; i < num_elements; i++) {
			PrintValue(h_data[i]);
			printf(", ");
		}
		printf("\n\n");
	}

	// Check
	int retval = CompareResults(h_data, h_reference, num_elements, verbose);

	// Cleanup
	if (h_reference) free(h_reference);
	if (h_data) free(h_data);

	return retval;
}


/**
 * Verify the contents of a device array match those
 * of a host array
 */
template <typename T>
void DisplayDeviceResults(
	T *d_data,
	size_t num_elements)
{
	// Allocate array on host
	T *h_data = (T*) malloc(num_elements * sizeof(T));

	// Reduction data back
	cudaMemcpy(h_data, d_data, sizeof(T) * num_elements, cudaMemcpyDeviceToHost);

	// Display data
	printf("\n\nData:\n");
	for (int i = 0; i < num_elements; i++) {
		PrintValue(h_data[i]);
		printf(", ");
	}
	printf("\n\n");

	// Cleanup
	if (h_data) free(h_data);
}



/******************************************************************************
 * Timing
 ******************************************************************************/


struct CpuTimer
{
#if defined(_WIN32) || defined(_WIN64)

	LARGE_INTEGER ll_freq;
	LARGE_INTEGER ll_start;
	LARGE_INTEGER ll_stop;

	CpuTimer()
	{
		QueryPerformanceFrequency(&ll_freq);
	}

	void Start()
	{
		QueryPerformanceCounter(&ll_start);
	}

	void Stop()
	{
		QueryPerformanceCounter(&ll_stop);
	}

	float ElapsedMillis()
	{
		double start = double(ll_start.QuadPart) / double(ll_freq.QuadPart);
		double stop  = double(ll_stop.QuadPart) / double(ll_freq.QuadPart);

		return (stop - start) * 1000;
	}

#else

	rusage start;
	rusage stop;

	void Start()
	{
		getrusage(RUSAGE_SELF, &start);
	}

	void Stop()
	{
		getrusage(RUSAGE_SELF, &stop);
	}

	float ElapsedMillis()
	{
		float sec = stop.ru_utime.tv_sec - start.ru_utime.tv_sec;
		float usec = stop.ru_utime.tv_usec - start.ru_utime.tv_usec;

		return (sec * 1000) + (usec / 1000);
	}

#endif
};

struct GpuTimer
{
	cudaEvent_t start;
	cudaEvent_t stop;

	GpuTimer()
	{
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
	}

	~GpuTimer()
	{
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
	}

	void Start()
	{
		cudaEventRecord(start, 0);
	}

	void Stop()
	{
		cudaEventRecord(stop, 0);
	}

	float ElapsedMillis()
	{
		float elapsed;
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&elapsed, start, stop);
		return elapsed;
	}
};



}// namespace b40c
