#include<iostream>
#include <fstream>
#include <string.h>
#include "file_read.h"
using namespace std;


int main(int argc, char **argv) {
    string file_name;

    int i = 1;
    while(i < argc) {
        if (strcmp( argv[i], "-F") == 0){
            file_name = argv[i+1];
        }
        i++;
    }


    vector<string> sequences = seqs_read(file_name);
    for (string i: sequences) {
        cout << i << endl;
    }


}
