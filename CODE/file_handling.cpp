/********** FUNCTIONS FOR FILE HANDLING ********/

#include "definitions.h"

using namespace std;

void Remove_last_line(string filename){

    /**** fstream object are used to read the lines in the file ****/
    fstream read_file;

    /**** Open the file with the provided filename ****/
    read_file.open(filename);

    /**** If file failed to open, exit with an error message and error exit status ****/
    if (read_file.fail()){
        cout << "Error opening file." << endl; 
    }

    /**** Create a vector to store all the file lines, and a string line to store each line that we read ****/ 
    vector<string> lines;
    string line;
    while (getline(read_file, line)){
        lines.push_back(line);
    }

    /**** The vector will now contain an element for each line in the file ****/
    read_file.close();

    /**** Create ofstream object for writing to the file */
    ofstream write_file;
    write_file.open(filename);

    /**** Write all the lines except the last one ****/
    for (int i = 0; i < lines.size() - 1; i++){
        write_file << lines[i] << endl;
    }

    /**** Close the filled file ****/
    write_file.close();
}