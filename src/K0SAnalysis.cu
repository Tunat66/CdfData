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
{}

// Start method
void K0SAnalysis::start() {
    // Currently empty, can be filled if needed
}

// Stop method
void K0SAnalysis::stop() {
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
    double* cPrimaryVertex = new double[2];
    cPrimaryVertex[0] = primaryVertex[0];
    cPrimaryVertex[1] = primaryVertex[1];
    //allocate memory on the GPU for the primary vertex
    double* d_primaryVertex;
    cudaMalloc(&d_primaryVertex, 2 * sizeof(double));
    // Copy the primary vertex to the GPU
    cudaMemcpy(d_primaryVertex, cPrimaryVertex, 2 * sizeof(double), cudaMemcpyHostToDevice);


    // copy over the cutoff values to the GPU
    // Set the cutoff values in constant memory
    //CUDA_CHECK(cudaMemcpyToSymbol(pTcutoff, &this->pTcutoff, sizeof(double)));
    //CUDA_CHECK(cudaMemcpyToSymbol(TrackImpactParametercutoff, &this->TrackImpactParametercutoff, sizeof(double)));
    //CUDA_CHECK(cudaMemcpyToSymbol(Lxycutoff, &this->Lxycutoff, sizeof(double)));
    //CUDA_CHECK(cudaMemcpyToSymbol(ImpactParametercutoff, &this->ImpactParametercutoff, sizeof(double)));
    
    // Allocate memory on the GPU
    double* d_massArray;
    double* d_lifetimeArray;
    int* d_massCounter;
    cudaMalloc(&d_massArray, numTracks * sizeof(double));
    cudaMalloc(&d_lifetimeArray, numTracks * sizeof(double));
    cudaMalloc(&d_massCounter, sizeof(int));

    // Copy track data to the GPU
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

    // Allocate memory for the flattened data on the device
    double* d_flattenedTrackData;
    cudaMalloc(&d_flattenedTrackData, flattenedTrackData.size() * sizeof(double));
    cudaMemcpy(d_flattenedTrackData, flattenedTrackData.data(), flattenedTrackData.size() * sizeof(double), cudaMemcpyHostToDevice);

    // Allocate memory for the array of pointers on the device
    const double** d_trackData;
    cudaMalloc(&d_trackData, numTracks * sizeof(double*));

    // Create an array of pointers on the host
    std::vector<double*> h_trackData(numTracks);
    size_t offset = 0;
    for (int i = 0; i < numTracks; ++i) {
        h_trackData[i] = d_flattenedTrackData + offset;
        offset += trackData[i].size();
    }

    // Copy the array of pointers to the device
    cudaMemcpy(d_trackData, h_trackData.data(), numTracks * sizeof(double*), cudaMemcpyHostToDevice);

    // Initialize massCounter on the GPU
    int massCounter = 0;
    cudaMemcpy(d_massCounter, &massCounter, sizeof(int), cudaMemcpyHostToDevice);
    // Launch the kernel
    dim3 blockSize(16, 16);  // 16x16 threads per block
    dim3 gridSize((numTracks + blockSize.x - 1) / blockSize.x,
              (numTracks + blockSize.y - 1) / blockSize.y);
    processTracksKernel<<<gridSize, blockSize>>>(d_trackData, numTracks, d_primaryVertex, d_massArray, d_lifetimeArray, d_massCounter, m_pion);
    // Check for kernel launch errors
    cudaGetLastError();
    // Synchronize the device to ensure the kernel has finished
    cudaDeviceSynchronize();


    // Copy results back to the host
    std::vector<double> massArray(numTracks);
    
    std::vector<double> lifetimeArray(numTracks);
    cudaMemcpy(massArray.data(), d_massArray, numTracks * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(lifetimeArray.data(), d_lifetimeArray, numTracks * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&massCounter, d_massCounter, sizeof(int), cudaMemcpyDeviceToHost);

    // Free GPU memory
    cudaFree(d_trackData);
    cudaFree(d_flattenedTrackData);
    cudaFree(d_massArray);
    cudaFree(d_lifetimeArray);
    cudaFree(d_massCounter);
    cudaFree(d_primaryVertex);

    // Free CPU memory
    delete[] cPrimaryVertex;
    // No need to free d_trackData with delete[] as it was allocated using cudaMalloc.
    // The cudaFree call above already handles the deallocation of d_trackData.

    // Store results in the class members
    this->massArray.insert(this->massArray.end(), massArray.begin(), massArray.end());
    this->lifetimeArray.insert(this->lifetimeArray.end(), lifetimeArray.begin(), lifetimeArray.end());
    //std::cout << "Mass counter: " << massCounter << std::endl;
    this->massCounter += massCounter;

    
}
}  // namespace Cdf

