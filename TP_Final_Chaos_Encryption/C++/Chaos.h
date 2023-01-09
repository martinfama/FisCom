//An object that uses Transmissor.h and Receptor.h to encrypt and
//decrypt WAV files

#ifndef CHAOS_H
#define CHAOS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include "Transmissor.h"
#include "Receptor.h"
#include "AudioFile.h" //library taken from: https://github.com/adamstark/AudioFile/

class Chaoscrypt {

    public:
        //this initialization opens an existing .chaos file and imports its data
        Chaoscrypt(std::string filename_) {
            filename = filename_;
            file.open(filename);
            if (!file.is_open()) {
                std::cout << "Error opening file: " << filename << "\n";
                return;
            }
            std::cout << "Opened file: " << filename << "\n" << std::flush;
            std::cout << "Reading data from file... " << std::flush;
            file >> t_grace;
            file >> t_tot;
            file >> freqSamp;
            file >> sampleRatio;
            file >> epsilon;
            double tmp_s;
            while (file >> tmp_s) {
                s.push_back(tmp_s);
            }
            std::cout << "Done\n" << std::flush;
        }
        //this initialization generates the option to make a new .chaos files
        Chaoscrypt(double t_g_, double freqSamp_, int sampleRatio_, double epsilon_) {
            t_grace = t_g_;
            t_tot = t_grace;
            freqSamp = freqSamp_;
            sampleRatio = sampleRatio_;
            epsilon = epsilon_;
        }

        void printData() {
            std::cout << "T_g:         " << t_grace << "\n";
            std::cout << "T_tot:       " << t_tot << "\n";
            std::cout << "FreqSamp:    " << freqSamp << "\n";
            std::cout << "SampleRatio: " << sampleRatio << "\n";
            std::cout << "Epsilon:     " << epsilon << "\n\n";
        }

        //output encrypted signal s(t) as wavfile
        void outputWAV(std::string wavFilename) {
            std::cout << "Outputting encrypted data to wavfile: " << wavFilename << "\n" << std::flush;
            AudioFile<double> audiofile;
            audiofile.setNumChannels(1);
            audiofile.setBitDepth(16);
            audiofile.setSampleRate(freqSamp);

            std::cout << "Manipulating wavData... " << std::flush;
            
            std::vector<double> tmp_wavData = s;
            tmp_wavData.erase(tmp_wavData.begin(), tmp_wavData.begin() + t_grace*freqSamp);
            
            //get samples every sampleRatio step
            std::vector<double> wavData;
            for (int i = 0; i < tmp_wavData.size(); i += sampleRatio) {
                wavData.push_back(tmp_wavData[i]);
            }
            
            //normalize signal 
            double max_wavData = std::max(*std::max_element(std::begin(wavData), std::end(wavData)), 
                                     std::abs(*std::min_element(std::begin(wavData), std::end(wavData))));
            for (int i = 0; i < wavData.size(); i++) {
                wavData[i] /= max_wavData;
            }
            std::cout << "Done\n" << std::flush;
            audiofile.samples[0] = wavData;
            std::cout << "Saving to file... " << std::flush; 
            audiofile.save(wavFilename);
            std::cout << "Done\n" << std::flush;
        }

        //receives wavFilename as an input WAV file, encrypts it and saves to a .chaos file
        void encryptWAV(std::string chaosFilename, std::string wavFilename) {
            std::cout << "Encrypting wavfile: " << wavFilename << " ...\n" << std::flush;
            AudioFile<double> wavData;
            wavData.load(wavFilename);
            std::cout << "Starting Lorenz system simulation... " << std::flush;
            transmissor t(t_grace, wavData.getSampleRate(), wavData.samples[0], sampleRatio, epsilon);
            s = t.run();
            std::cout << "Done\n" << std::flush;
            t_tot = s.size()/freqSamp;
            std::cout << "Writing to chaosfile: " << chaosFilename << " ... " << std::flush;
            std::ofstream chaosFile(chaosFilename);
            chaosFile << t_grace << "\n";
            chaosFile << t_tot << "\n";
            chaosFile << freqSamp << "\n";
            chaosFile << sampleRatio << "\n";
            chaosFile << epsilon << "\n";
            for (int i = 0; i < s.size(); i++) {
                chaosFile << s[i] << " ";
            }
            chaosFile.close();
            std::cout << "Done\n\n";
        }

        //decrypts signal s(t)
        void decryptToWAV(std::string wavFilename) {
            std::cout << "Decrypting self...\n" << std::flush;
            AudioFile<double> audiofile;
            receptor r(freqSamp, s, epsilon);
            audiofile.setNumChannels(1);
            audiofile.setBitDepth(16);
            audiofile.setSampleRate(freqSamp);
            std::cout << "Starting Lorenz simulation... " << std::flush;
            std::vector<double> tmp_wavData = r.run();
            std::cout << "Done\n" << std::flush;
            std::cout << "Manipulating wavData... " << std::flush;
            tmp_wavData.erase(tmp_wavData.begin(), tmp_wavData.begin() + t_grace*freqSamp);
            std::vector<double> wavData;
            for (int i = 0; i < tmp_wavData.size(); i += sampleRatio) {
                wavData.push_back(tmp_wavData[i]);
            }
            
            //normalize. but only if data exceeds -1 or 1
            double max_wavData = std::max(*std::max_element(std::begin(wavData), std::end(wavData)), 
                                     std::abs(*std::min_element(std::begin(wavData), std::end(wavData))));
            if (max_wavData > 1) { 
                for (int i = 0; i < wavData.size(); i++) {
                    wavData[i] /= max_wavData;
                }
            }
            std::cout << "Done\n" << std::flush;
            audiofile.samples[0] = wavData;
            std::cout << "Saving to wavfile: " << wavFilename << " ... " << std::flush;
            audiofile.save(wavFilename);
            std::cout << "Done\n" << std::flush;
        }

    private:
        double t_grace; //time given to systems to synchronize
        double t_tot;   //total time of signal
        double freqSamp = 44100; //sampling frequency
        int sampleRatio; //atractor sampling distance
        double epsilon; //message amplitude modulation
        std::vector<double> s; //signal
        std::ifstream file;
        std::string filename;
};

#endif