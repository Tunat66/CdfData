#include "K0SAnalysis.hh"
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << " - " \
                      << cudaGetErrorString(err) << std::endl; \
            exit(EXIT_FAILURE); \
        } \
    } while (0)
namespace Cdf 
{
// Constructor
K0SAnalysis::K0SAnalysis(double pTcutoff, double TrackImpactParametercutoff, double Lxycutoff, double ImpactParametercutoff)
    : pTcutoff(pTcutoff),
      TrackImpactParametercutoff(TrackImpactParametercutoff),
      Lxycutoff(Lxycutoff),
      ImpactParametercutoff(ImpactParametercutoff) 
{
    cPrimaryVertex = new double[2];
}

// Start method
void K0SAnalysis::start() {
    // Currently empty, can be filled if needed
    // Allocate GPU memory once
    cudaMalloc(&d_massArray, maxPairs * sizeof(double));
    cudaMalloc(&d_lifetimeArray, maxPairs * sizeof(double));
    
    cudaMalloc(&d_massCounter, sizeof(int));
    cudaMalloc(&d_primaryVertex, 2 * sizeof(double));

    cudaMalloc(&d_flattenedTrackData, maxTracksPerEvent * 5 * sizeof(double));  // Assuming 5 parameters per track
    cudaMalloc(&d_trackData, maxTracksPerEvent * sizeof(double*));
}

// Stop method
void K0SAnalysis::stop() {

    
    
    //free CPU memory
    delete[] cPrimaryVertex;


    // Free GPU memory
    cudaFree(d_trackData);
    cudaFree(d_flattenedTrackData);
    cudaFree(d_massArray);
    cudaFree(d_lifetimeArray);
    cudaFree(d_massCounter);
    cudaFree(d_primaryVertex);

    // Remove zeroes from the mass and lifetime arrays
    massArray.erase(std::remove(massArray.begin(), massArray.end(), 0), massArray.end());
    lifetimeArray.erase(std::remove(lifetimeArray.begin(), lifetimeArray.end(), 0), lifetimeArray.end());

    //put the data out into Dat files
    // Write massArray to a file
    // Remove zeroes and NaN entries from the mass and lifetime arrays
    massArray.erase(std::remove_if(massArray.begin(), massArray.end(),
                                   [](double value) { return value == 0 || std::isnan(value); }),
                    massArray.end());
    lifetimeArray.erase(std::remove_if(lifetimeArray.begin(), lifetimeArray.end(),
                                       [](double value) { return value == 0 || std::isnan(value); }),
                        lifetimeArray.end());
    std::ofstream massFile("massArray.dat");
    if (massFile.is_open()) {
        for (const auto& mass : massArray) {
            massFile << mass << "\n";
        }
        massFile.close();
        std::cout << "Mass data written to massArray.dat" << std::endl;
    } else {
        std::cerr << "Error: Could not open massArray.dat for writing." << std::endl;
    }

    // Write lifetimeArray to a file
    std::ofstream lifetimeFile("lifetimeArray.dat");
    if (lifetimeFile.is_open()) {
        for (const auto& lifetime : lifetimeArray) {
            lifetimeFile << lifetime << "\n";
        }
        lifetimeFile.close();
        std::cout << "Lifetime data written to lifetimeArray.dat" << std::endl;
    } else {
        std::cerr << "Error: Could not open lifetimeArray.dat for writing." << std::endl;
    }

    

}

// Event method
void K0SAnalysis::event(Event* ev) {
    // Get the number of tracks
    auto tracks = ev->getTracks();
    int numTracks = tracks.size();
    std::array<double, 2> primaryVertex = ev->getVertex();
    
    cPrimaryVertex[0] = primaryVertex[0];
    cPrimaryVertex[1] = primaryVertex[1];
    //allocate memory on the GPU for the primary vertex
    // Copy the primary vertex to the GPU
    cudaMemcpy(d_primaryVertex, cPrimaryVertex, 2 * sizeof(double), cudaMemcpyHostToDevice);


    // copy over the cutoff values to the GPU
    // Set the cutoff values in constant memory
    //CUDA_CHECK(cudaMemcpyToSymbol(pTcutoff, &this->pTcutoff, sizeof(double)));
    //CUDA_CHECK(cudaMemcpyToSymbol(TrackImpactParametercutoff, &this->TrackImpactParametercutoff, sizeof(double)));
    //CUDA_CHECK(cudaMemcpyToSymbol(Lxycutoff, &this->Lxycutoff, sizeof(double)));
    //CUDA_CHECK(cudaMemcpyToSymbol(ImpactParametercutoff, &this->ImpactParametercutoff, sizeof(double)));
    

    // Gather Track Data using the event object
    std::vector<std::vector<double>> trackData(numTracks);
    for (int i = 0; i < numTracks; ++i) {
        trackData[i] = tracks[i].getTrackParameters(); // Replace with actual track parameter
    }
    //each track has 5 pieces of data, so we need to allocate 5 times the number of tracks
    //see Track.hh for the pieces of the data (protected attributes)

    // Flatten the track data into a single array
    std::vector<double> flattenedTrackData;
    for (const auto& track : trackData) {
        flattenedTrackData.insert(flattenedTrackData.end(), track.begin(), track.end());
    }

    // we will keep reusing the d_flattenedTrackData and d_trackData pointers
    // no need to clear it as it will be overwritten and numTracks will be passed
    // to the kernel to prevent it from reading junk data

    // Allocate memory for the flattened data on the device
    cudaMemcpy(d_flattenedTrackData, flattenedTrackData.data(), flattenedTrackData.size() * sizeof(double), cudaMemcpyHostToDevice);
    // Create an array of pointers on the host
    if (numTracks > maxTracksPerEvent) {
        std::cerr << "Error K0SAnalysis: numTracks exceeds maxTracksPerEvent!" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::array<double*, maxTracksPerEvent> h_trackData; 
    size_t offset = 0;
    for (int i = 0; i < numTracks; ++i) {
        h_trackData[i] = d_flattenedTrackData + offset; //doing address arithmetic
        offset += trackData[i].size();
    }
    // Copy the array of pointers to the device
    cudaMemcpy(d_trackData, h_trackData.data(), numTracks * sizeof(double*), cudaMemcpyHostToDevice);

    // Launch the kernel
    cudaMemset(d_massCounter, 0, sizeof(int));
    cudaMemset(d_massArray, -1, maxPairs * sizeof(double)); //set masses to invalid values
    dim3 blockSize(16, 16);  // 16x16 threads per block
    dim3 gridSize((numTracks + blockSize.x - 1) / blockSize.x,
              (numTracks + blockSize.y - 1) / blockSize.y);
    processTracksKernel<<<gridSize, blockSize>>>(d_trackData, numTracks, d_primaryVertex, d_massArray, d_lifetimeArray, d_massCounter, m_pion);
    // Check for kernel launch errors
    cudaGetLastError();
    // Synchronize the device to ensure the kernel has finished
    cudaDeviceSynchronize();

    //now extract the data from d_massArray
    // copy the amount of masses that were calculated
    int massCounter = 0;
    cudaMemcpy(&massCounter, d_massCounter, sizeof(int), cudaMemcpyDeviceToHost);
    //std::cout << "Mass counter: " << massCounter << std::endl;
    if (massCounter > maxPairs) {
        std::cerr << "Error: massCounter exceeds allocated size of d_massArray!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Append the values in d_massArray to the massArray
    std::vector<double> tempMassArray(maxPairs);
    std::vector<double> tempLifetimeArray(maxPairs);

    

    cudaMemcpy(tempMassArray.data(), d_massArray, maxPairs * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(tempLifetimeArray.data(), d_lifetimeArray, maxPairs * sizeof(double), cudaMemcpyDeviceToHost);

    //clean tempMassArray and tempLifetimeArray of invalid values
    tempMassArray.erase(std::remove_if(tempMassArray.begin(), tempMassArray.end(),
                                       [](double value) { return value == -1 || std::isnan(value); }),
                        tempMassArray.end());

    //for (int i = 0; i < massCounter; ++i) {
    //    std::cout << "Mass: " << tempMassArray[i] << std::endl;
    //}
    massArray.insert(massArray.end(), tempMassArray.begin(), tempMassArray.end());
    lifetimeArray.insert(lifetimeArray.end(), tempLifetimeArray.begin(), tempLifetimeArray.end());

}
}  // namespace Cdf

