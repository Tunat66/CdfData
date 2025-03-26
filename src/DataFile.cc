#include "DataFile.hh"

namespace Cdf
{
    DataFile::DataFile(std::string filename)
    {
        //strong typing allows me to assign straight away
        //OPEN THE FILE AND READ
        fileStream.open(filename, std::ios::in);
        if (!fileStream.is_open())
        {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            throw std::ios_base::failure("File open failed");
        }
    }
    DataFile::~DataFile()
    {
        if (fileStream.is_open())
        {
            fileStream.close();
        }
    }


    //sequential access to prevent breaking the system
    Event* DataFile::next()
    {
        std::string line;
        Event* event = new Event();
        // Read the next valid line from the file
        // while loop is needed to skip any comments
        while (std::getline(fileStream, line))
        {
            //trim the line and apply skip logic
            line = trim(line);
            if (!skipLogic(line))
            {
                continue;
            }
            
            //else parse the event header
            std::istringstream headerStream(line);
            int runNumber, eventNumber, numTracks;
            double vertexX, vertexY;
            //std::cout << "Parsing event header: " << line << std::endl;
            if (!(headerStream >> runNumber >> eventNumber >> vertexX >> vertexY >> numTracks)) {
                std::cerr << "Error: Failed to parse event header: " << line << std::endl;
                return nullptr; // Unexpected header
            }
            // otherwise read the tracks in the event
            std::vector<Track> tracks; // Create a vector to store Track objects
            std::vector<std::vector<double>> trackData(numTracks, std::vector<double>(5));
            for (int i = 0; i < numTracks; ++i) {
                double p1, p2, p3, p4, p5; // track parameters
                //std::cout << "Parsing track " << i + 1 << " of " << numTracks << std::endl;
                if (!(fileStream >> p1 >> p2 >> p3 >> p4 >> p5)) {return nullptr;}
                // Create a Track object and push it into the vector
                std::vector <double> trackParameters = { p1, p2, p3, p4, p5 };
                Track track(trackParameters);
                tracks.emplace_back(track);
            }
            // populate the event
            event->runNumber = runNumber;
            event->eventNumber = eventNumber;
            event->vertex[0] = vertexX;
            event->vertex[1] = vertexY;
            event->tracks = tracks;
            // Consume the remaining newline character after track data
            fileStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            return event;
        }
        return nullptr;
    }
    //create event bundles to be processed, ideally 256 events long
    std::vector<Event*> DataFile::next_eventBundle()
    {
        std::vector<Event*> bundle = {};
        for (int i = 0; i < BundleSize; ++i)
        {
            bundle.push_back(next());
        }
        return bundle;
    }
    //these event bundles are then passed to the gpu kernels to do stuff

    
    //UTILITY STUFF
    // Rewind the file to the beginning
    void DataFile::rewind() {
        if (fileStream.is_open()) {
            fileStream.clear(); // Clear EOF or error flags
            fileStream.seekg(0, std::ios::beg); // Move the file pointer to the beginning
        }
    }

    std::string DataFile::trim(std::string line)
    {
        // Trim leading/trailing whitespace
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);
        return line;
    }

    bool DataFile::skipLogic(std::string line)
    {
        //add any skip criteria here
        if (line.empty()) {
            return false; //skip empty lines
            std::cout << "Warning: Empty line in .dat file!" << line << std::endl;
        }
        if (line[0] == '#') {
            //skip comment lines
            std::cout << line << std::endl;
            return false; 
        }
        //return true otherwise
        return true;
    }    

}