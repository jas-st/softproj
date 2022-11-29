#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <sstream>
using namespace std;

vector<string> seqs_read(string file_name) {
    vector<string> sequences;
    string line;
    string seq;
    bool first_line = true;

    ifstream seq_file(file_name);
    if (!seq_file) {
        cout <<"Unable to open file";
        exit(1);
    }

    while(getline(seq_file, line)) {
        if (first_line) {
            first_line = false;
        }
        else {
            istringstream s(line);
            for (int i = 0; i<2; i++) {
                s >> seq;
            }
            sequences.push_back(seq);
        }
    }

    return sequences;

}
