/*
   Function: remove 16 bp barcode and 6 bp linker sequence from read 1.
   Written by Hequan Sun, MPIPZ, Email:sunhequan@gmail.com/sun@mpipz.mpg.de
*/
#include <iostream>
#include  <fstream>
#include  <sstream>
#include      <map>
#include   <vector>
#include   <string>
#include <stdlib.h>
#include <algorithm>
#include <sys/stat.h>
#include <dirent.h>
#include <assert.h>
#include <time.h>  /* clock_t, clock, CLOCKS_PER_SEC */
#include "split_string.h"
// zlib library for ungzip/gzip in c++ style
#include   "./gzlib/gzstream.h"
// ogzstream
// igzstream
bool trim_reads(char* fastqfilename1, 
                    char* fastqfilename2, 
                    int   barcode_len);
bool create_folder(string tmpfolder);
using namespace std;
int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        cout << "\nFunction: trim 16bp barcode + 6 linker seqs off R1 reads of paired-end fastq.gz files -- ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\nUsage: 10X_barcode_trimmer this_R1.fastq.gz this_R2.fastq.gz"  << endl << endl;
        cout << "      barcode_len by default is 16; linker length 6 bp. So 22 bp off R1 reads. " << endl;
        cout << "      this_R1.fastq.gz this_R2.fastq.gz: read files -- only gzipped .gz files accepted" << endl << endl;
        cout << "Note1: pls keep order of inputs as given"                        << endl;
        cout << "Note2: output will be trimmed_R1.fastq.gz trimmed_R2.fastq.gz."  << endl << endl;
        return 1;
    }
    clock_t tbeg;
    tbeg = clock();
    
    int barcode_len  = 16;
    if(!trim_reads(argv[1], argv[2], barcode_len) )
    {
        return false;
    }
    cout << "\n  Total time consumed : " << (float)(clock()-tbeg)/CLOCKS_PER_SEC << " seconds.\n" << endl;
    return 0;
}

bool trim_reads(char* fastqfilename1, char* fastqfilename2, int barcode_len)
{
    clock_t tbeg;
    tbeg = clock();
    // first in file
    igzstream fp1(fastqfilename1);
    if(!fp1.good())
    {
        cout << "Cannot open file "     << fastqfilename1 << " to read reads. Exited.\n";
        return false;
    }
    else
    {
        cout << "   Info: R1 file: "    << fastqfilename1 << endl;
    }
    // first output file
    string ss1;
    vector<string> fastqinfo1 = split_string(fastqfilename1, '/');
    ss1 = "trimmed_" + fastqinfo1[fastqinfo1.size()-1];    
    ogzstream ofp1(ss1.c_str());
    // second in file
    igzstream fp2(fastqfilename2);
    if(!fp2.good())
    {
        cout << "Cannot open file "     << fastqfilename2 << " to read reads. Exited.\n";
        return false;
    }
    else
    {
        cout << "   Info: R2 file: "    << fastqfilename2 << endl;
    }    
    // second output file
    string ss2;
    vector<string> fastqinfo2 = split_string(fastqfilename2, '/');
    ss2 = "trimmed_" + fastqinfo2[fastqinfo2.size()-1];
    ogzstream ofp2(ss2.c_str());    
    //
    cout << "\n  Trimming 16 bp barcode and 6 bp linker off R1 reads... " << endl;        
    // parse the original fastq.gz file
    unsigned long checki   = 0;
    while(fp1.good() && fp2.good())
    {
        string line1("");
        getline(fp1, line1);
        if(line1.size()==0) continue;
        //
        string line2("");
        getline(fp2, line2);
        if(line2.size()==0) continue;
        
        if(line1.substr(0, 1).compare("@")==0) // first line with @XXX is the read id
        {  
            // output current read
            checki ++;
            if(checki%50000 == 0) 
            cout << "   Info: " << checki << " read pairs checked.." << endl;
            string seq1(""); getline(fp1, seq1);
            string plu1(""); getline(fp1, plu1); 
            string qua1(""); getline(fp1, qua1);            
            string this_barcode = seq1.substr(0, barcode_len);
            //
            string seq2(""); getline(fp2, seq2);
            string plu2(""); getline(fp2, plu2);
            string qua2(""); getline(fp2, qua2);
            // out R1: note: remove 16 bp plus 6 bp
            ofp1 << line1 << " bc:"  << this_barcode << " UMI:" << seq1.substr(16, 6) << endl
                 << seq1.substr(22)  << endl
                 << plu1             << endl
                 << qua1.substr(22)  << endl;
            // out R2
            ofp2 << line2 << " bc:"  << this_barcode << endl
                 << seq2  << endl
                 << plu2  << endl
                 << qua2  << endl;
        }
    }
    cout << "   Info: " << checki << " read pairs trimmed in total." << endl;    
    fp1.close();
    fp2.close();
    ofp1.close();
    ofp2.close();
    //
    cout << "  Time consumed in sampling : " << (float)(clock()-tbeg)/CLOCKS_PER_SEC << " seconds." << endl;
    return true;
}
//
bool create_folder(string tmpfolder)
{
    DIR* dir = opendir(tmpfolder.c_str());
    if (dir)
    {
        // Directory exists.
        closedir(dir);                
    }
    else if (ENOENT == errno)
    {
        // Directory does not exist.
        const int dir_err = mkdir(tmpfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (dir_err == -1)
        {
            cout << "   Error: cannot create directory " << tmpfolder << endl;
            return false;
        }
    }
    else;
    //
    return true;
}

