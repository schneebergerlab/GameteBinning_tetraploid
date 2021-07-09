/* 
this function converts a readname-sorted 10x bam (from cellranger-dna cnv outpt) to fastq, with corrected barcodes (raw barcodes in readname).
*/
#include   <iostream>
#include    <fstream>
#include    <sstream>
#include        <map>
#include     <vector>
#include     <string>
#include   <stdlib.h>
#include  <algorithm>
#include <sys/stat.h>
#include   <dirent.h>
#include   <assert.h>
#include     <time.h>  /* clock_t, clock, CLOCKS_PER_SEC */
#include "split_string.h"
#include "./gzlib/gzstream.h"
int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        cout << "\n Function: convert 10x bam from \"cellranger-dna cnv\" to paired-end fq.gz:\n"
             << " Usage: samtools view RNsorted_bam.bam | T10xbam2fq - out_prefix_str " << endl << endl;  
        return false;
    }
    clock_t tbeg;
    tbeg = clock();
    string outprefix = (string)argv[2];
    if(outprefix.size()==0) outprefix = "fun_convert";
    //
    bool verbose = true; // to add future option
    // open files under folder
    stringstream convertedfqfile1;
    convertedfqfile1.str("");
    convertedfqfile1 << outprefix << "_fqfrom10xBam_bxCorrected_R1.fq.gz";
    ogzstream of1;
    of1.open((convertedfqfile1.str()).c_str());   
    if(!of1.good())
    {
         cout << "Cannot open new file " << convertedfqfile1.str() << " to write R1 reads. Exited.\n";
         return false;
    }
    //
    stringstream convertedfqfile2;
    convertedfqfile2.str("");
    convertedfqfile2 << outprefix << "_fqfrom10xBam_bxCorrected_R2.fq.gz";
    ogzstream of2;
    of2.open((convertedfqfile2.str()).c_str());   
    if(!of2.good())
    {
         cout << "Cannot open new file " << convertedfqfile2.str() << " to write R2 reads. Exited.\n";
         return false;
    }
    //    
    cout << "   Info: reads with corrected barcode will be collected in " << endl
         << "      " << convertedfqfile1.str() << endl
         << "      " << convertedfqfile2.str() << endl;
    cout << "   Info: R1 readname-line would be followed by \"1:N:0:0\" and the raw barcode from sequencing; " 
         << " while the first 16bp in read 1 sequence is the cellranger-dna corrected barcode according to whitelist. " 
         << endl;
    //
    unsigned long numraw   = 0;
    unsigned long total_R1 = 0;
    unsigned long total_R2 = 0;
    map<string, int> paired_name;
    //
    std::string line;
    while (std::getline(std::cin, line)) 
    {
        if(line.compare("quit")==0 || line.compare("q")==0 || line.compare("exit")==0) break;
        if(line.size()==0 || line[0]=='#') continue;
        numraw ++;
        if(numraw%100000000 == 0)
        {
            cout << "   info: " << numraw << "th read..." << endl;
        }
        //
        vector<string> lineinfo = split_string(line, '\t');
        string readname = lineinfo[0];
        string readseq  = lineinfo[9];
        string basequa  = lineinfo[10];
        //
        /* CB:Z	Chromium cellular barcode sequence that is error-corrected and 
                confirmed against a list of known-good barcode sequences.
                The cell barcode CB tag includes a suffix "-1" 
                    that labels the GEMs from a single channel and we call a GEM group.
           CR:Z	Chromium cellular barcode sequence as reported by the sequencer.
           CY:Z	Chromium cellular barcode read quality. Phred scores as reported by sequencer.
           //
           The one after "CB:Z:" will be catted in the beginning of R1.
        */
        // get raw barcode with after info: "\tCR:Z:GATCTAGGTGCGGTAA"
        size_t pos1  = line.find("\tCR:Z:");
        if(pos1 == std::string::npos)
        {
            cout << "   Error: not CR:Z field found in line: " << line << endl;
            return 1;
        }
        pos1 += 6;
        size_t pos2  = line.find("\t", pos1);
        if(pos2 == std::string::npos)
        {
            pos2 = line.size();
        }        
        string rawbcod = line.substr(pos1, pos2-pos1);
        // get raw barcode base quality
        pos1  = line.find("\tCY:Z:");
        if(pos1 == std::string::npos)
        {
            cout << "   Error: not CY:Z: field found in line: " << line << endl;
            return 1;
        }
        pos1 += 6;
        pos2  = line.find("\t", pos1);
        if(pos2 == std::string::npos)
        {
            pos2 = line.size();
        }        
        string rawbcqua = line.substr(pos1, pos2-pos1);
        
        // get corrected barcode
        pos1  = line.find("\tCB:Z:");
        string corrbcod = rawbcod; // initialized as raw barcode
        if(pos1 != std::string::npos)
        {
            pos1 += 6; // get corrected barcode from 'CB:Z:barcode-1'
            pos2  = line.find("\t", pos1);
            if(pos2 == std::string::npos)
            {
                pos2 = line.size();
            }  
            corrbcod = line.substr(pos1, pos2-pos1-2); // excluding 2 additional positions for "-1"
        }
        //
        int hexflag = strtol(lineinfo[1].c_str(), NULL, 0);
        // get read as R1 or R2
        string readflag = "U";
        if((hexflag & 0x40) == 0x40 && (hexflag & 0x80) == 0)
        {
            readflag = "R1";
        }
        else
        if((hexflag & 0x40) == 0    && (hexflag & 0x80) == 0x80)
        {
            readflag = "R2";
        }
        else
        {
            readflag = "U0";
            cout << "   hexflag=" << hexflag << endl;
            cout << "   hexflag & 0x40 == 0x40= " << (hexflag & 0x40) << endl;
            cout << "   hexflag & 0x80 == 0   = " << (hexflag & 0x80) << endl;            
        }
        // get reverse or non-reverse state of the read
        string reversed = "nrc";
        if((hexflag & 0x10) == 0x10) reversed = "rc";
        if(reversed.compare("rc")==0)
        {
            string this_read_rc = readseq;                    
            /* reverse sequence    */
            std::reverse(this_read_rc.begin(), this_read_rc.end());
            /* lowercase sequence  */
            std::transform(this_read_rc.begin(), this_read_rc.end(),this_read_rc.begin(), ::tolower);
            /* complement sequence */
            std::replace(this_read_rc.begin(), this_read_rc.end(), 'a', 'T');
            std::replace(this_read_rc.begin(), this_read_rc.end(), 'c', 'G'); // caution1: here c:G
            std::replace(this_read_rc.begin(), this_read_rc.end(), 'g', 'C'); // caution2: here g:C
            std::replace(this_read_rc.begin(), this_read_rc.end(), 't', 'A'); // if caution 1&2 is not consistent,
            std::replace(this_read_rc.begin(), this_read_rc.end(), 'u', 'A'); // they will result in c:G+G:C=c:C
            //
            std::transform(this_read_rc.begin(), this_read_rc.end(),this_read_rc.begin(), ::toupper);
            //
            readseq = this_read_rc;
            // also reverse quality
            std::reverse(basequa.begin(), basequa.end());
        }
        // chr	pos	CIGAR	barcode	MI	read-ordering:R1/R2	reversed	read-id	pnext
        // caution!!! code below also works for reversed-cases of alignment!
        // in bam/sam, it is always pos + CIGAR"M" ==> spanning of a read!!!
        if(readflag.compare("R1")==0)
        {
            of1 << "@"       << readname << " 1:N:0:0 " << " CR:Z:" << rawbcod << " " << reversed << endl;
            of1 << corrbcod  << readseq  << endl
                << "+"                   << endl
                <<  rawbcqua << basequa  << endl;
            total_R1 ++;
            //
            if(paired_name.size()>0 && paired_name.find(readname) == paired_name.end()) // not paired with last
            {
                cout << "   Error: read " << readname << " is not paired with last " << (*(paired_name.begin())).first << endl;
                return 1;
            }else
            if(paired_name.size()>0 && paired_name.find(readname) != paired_name.end()) // paired with last
            {
                paired_name.clear();
            }else
            if(paired_name.size()==0)
            {
                paired_name.insert(std::pair<string, int>(readname, 1)); // new read
            }
        }else
        if(readflag.compare("R2")==0)
        {
            of2 << "@"      << readname << " 1:N:0:0 " << " CR:Z:" << rawbcod << " " << reversed << endl;
            of2 <<             readseq  << endl
                << "+"                  << endl
                <<             basequa  << endl;
            total_R2 ++;
            //
            if(paired_name.size()>0 && paired_name.find(readname) == paired_name.end()) // not paired with last
            {
                cout << "   Error: read " << readname << " is not paired with last " << (*(paired_name.begin())).first << endl;
                return 1;
            }else
            if(paired_name.size()>0 && paired_name.find(readname) != paired_name.end()) // paired with last
            {
                paired_name.clear();
            }else
            if(paired_name.size()==0)
            {
                paired_name.insert(std::pair<string, int>(readname, 1)); // new read
            }
        }else
        {
            cout << "   Error: reads neither R1 nor R2, weired?" << endl;
            return 1;
        }
        line.clear();
    }
    // close files
    of1.close();
    of2.close();
    // check info
    if(total_R1 != total_R2)
    {
        cout << "   Error: R1 and R2 have different numbers: " << "R1=" << total_R1 << " vs R2=" << total_R2 << endl;
        return 1;
    }
    cout << "   Info: in total " 
         << numraw  
         << " aligment lines, of which " 
         << total_R1               
         << " read pairs with barcode info coverted." 
         <<  endl;
    cout << "   Info: time on converting streaming bam info: " 
         << (float)(clock()-tbeg)/CLOCKS_PER_SEC 
         << " seconds.\n" 
         << endl;
    return 0;
}
