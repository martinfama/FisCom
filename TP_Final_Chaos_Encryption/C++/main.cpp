//Encrypt and decrypt WAV files using synchronization of Lorenz systems
//Allows user to encrypt a 16-bit mono WAV file. Encrypted file is of type .chaos
//The .chaos format is as follows:

//1  t_grace     = time given to systems to synchronize before adding message to signal
//2  t_tot       = total time of signal
//3  freqSamp    = sampling frequency of signal (and therefore WAV file. default = 44100)
//4  sampleRatio = whole positive number. indicates de sampling distance of atractor
//5  epsilon     = amplitude modulation of m(t). Sent signal will be s(t) = u(t) + epsilon*m(t) 
//6  s           = {s_0, s_1, s_2, ..., s_n} signal s(t_n), saved as a vector of doubles  

//compilation: g++ -std=c++11 -o main main.cpp
//make sure to keep AudioFile.h, Transmissor.h, Receptor.h and Chaos.h in same folder as main.cpp

#include <cstdlib>
#include <iostream>
#include <string>

#include "AudioFile.h"
#include "Transmissor.h"
#include "Receptor.h"
#include "Chaos.h"

using namespace std;

int main(int argc, char** argv) {

    if (argc == 1) {
        cout << "No arguments passed. Usage: \n\n";
        
        cout << "To encrypt a WAV file:\n";
        cout << "  ./chaos -encrypt -i <wavfile_in> -o <chaosfile_out> -s <N>\n";
        cout << "        wavfile_in = .wav to encrypt (must be 16-bit mono)\n";
        cout << "     chaosfile_out = .chaos file to write encrypted message\n";
        cout << "                 N = Whole number, sample distance\n\n";

        cout << "To decrypt a CHAOS file:\n";
        cout << "  ./chaos -decrypt -i <chaosfile_in> -o <wavfile_out>\n";
        cout << "      chaosfile_in = .chaos to decrypt\n";
        cout << "       wavfile_out = .wav file to write decrypted message\n\n";

        cout << "To write CHAOS file to WAV:\n";
        cout << "  ./chaos -outwav -i <chaosfile_in> -o <wavfile_out>\n";
        cout << "      chaosfile_in = .chaos to output\n";
        cout << "       wavfile_out = .wav file to write encrypted message\n\n";
        return 0;
    }

    string argv_str = argv[1];

    if (argv_str == "-encrypt") {
        string wavfile_in_name;
        string chaosfile_out_name;
        int sampleRatio;
        for (int i = 2; i < argc;) {
            argv_str = argv[i];
            if (argv_str == "-i") {
                wavfile_in_name = argv[i+1];
                i += 2;
            }
            else if (argv_str == "-o") {
                chaosfile_out_name = argv[i+1];
                i += 2;
            }
            else if (argv_str == "-s") {
                sampleRatio = atoi(argv[i+1]);
                i += 2;
            }
        }
        Chaoscrypt chaosfile(10, 44100, sampleRatio, 0.01);
        chaosfile.encryptWAV(chaosfile_out_name, wavfile_in_name);
        return 0;
    }

    if (argv_str == "-decrypt") {
        string wavfile_out_name;
        string chaosfile_in_name;
        for (int i = 2; i < argc;) {
            argv_str = argv[i];
            if (argv_str == "-i") {
                chaosfile_in_name = argv[i+1];
                i += 2;
            }
            else if (argv_str == "-o") {
                wavfile_out_name = argv[i+1];
                i += 2;
            }
        }
        Chaoscrypt chaosfile(chaosfile_in_name);
        chaosfile.decryptToWAV(wavfile_out_name);
        return 0;
    }

    if (argv_str == "-output") {
        string wavfile_out_name;
        string chaosfile_in_name;
        for (int i = 2; i < argc;) {
            argv_str = argv[i];
            if (argv_str == "-i") {
                chaosfile_in_name = argv[i+1];
                i += 2;
            }
            else if (argv_str == "-o") {
                wavfile_out_name = argv[i+1];
                i += 2;
            }
        }
        Chaoscrypt chaosfile(chaosfile_in_name);
        chaosfile.outputWAV(wavfile_out_name);
        return 0;
    }

    return 0;
}
