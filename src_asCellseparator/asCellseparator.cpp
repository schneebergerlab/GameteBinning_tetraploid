/*
   note: we remove the first 16 bp from Read 1 (and quality string)!
   update: remove 16 bp barcode and 6 bp linker sequence from read 1.
   Writen by Hequan Sun, MPIPZ
   Email: sun@mpipz.mpg.de/sunhequan@gmai.com
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
struct READ
{
    string R1;// name \n seq \n + \n quality
    string R2;
};
map<string, vector<READ> > myBCReads; // barcode separated reads
map<string, vector<READ> >::iterator myBCitr;
int  effective_cutoff = 500;          // caution: self-defined, need tuning for future
int  min_readpair     = 500000;       // caution: self-defined, need tuning for future
string output_folder  = "asCellseparator_sep_cells";
//
bool create_folder(string tmpfolder);
bool get_barcode_stat(char* fastqfilename1, int barcode_len, int effective_cutoff, map<string, unsigned long>* barcode_stat);
bool output_tmp_reads(map<string, vector<READ> >* myBCReads, int ipart); // only output those with 100 more reads
bool separate_reads(char* fastqfilename1, char* fastqfilename2, int barcode_len, map<string, unsigned long> barcode_stat);
using namespace std;
int main(int argc, char* argv[])
{
    if(argc < 7)
    {
        cout << "\nFunction: separate reads from paired-end fastq.gz files with self-detected barcodes -- ";
        const char *buildString = __DATE__ ", " __TIME__ ".";
        cout << "(compiled on " << buildString << ")" << endl;
        cout << "\nUsage: asCellseparator barcode_len this_R1.fastq.gz this_R2.fastq.gz effectiveReadPair min_readpair output_folder"  << endl << endl;
        cout << "      barcode_len                      : length of 10X barcodes to extract"                       << endl;
        cout << "      this_R1.fastq.gz this_R2.fastq.gz: read files -- only gzipped .gz files accepted"           << endl;
        cout << "      effectiveReadPair                : minimum number of read pairs to define a good barcode. " << endl << endl;
        cout << "      min_readpair                     : minimum number of read pairs to output as a part. "      << endl << endl;        
        cout << "Note1: pls keep order of inputs as given"            << endl;
        cout << "Note2: output will be partI_barcode_R1.fastq.gz partI_barcode_R2.fastq.gz." << endl << endl;
        return 1;
    }
    clock_t tbeg;
    tbeg = clock();    
    int barcode_len  = atoi(argv[1]);
    effective_cutoff = atoi(argv[4]);
    min_readpair     = atoi(argv[5]);
    output_folder    = (string)argv[6];
    if(!create_folder(output_folder))
    {
        return false;
    }
    map<string, unsigned long> barcode_stat; // only collect barcodes with read pairs >= effective_cutoff
    if(!get_barcode_stat(argv[2], barcode_len, effective_cutoff, &barcode_stat))
    {
        return 1;
    } 
    cout << "   Info: number of effective barcodes: " << barcode_stat.size() << endl;
    cout << "   Total time consumed for checking barcode stat: " 
         << (float)(clock()-tbeg)/CLOCKS_PER_SEC << " seconds.\n" << endl;       
    //
    if(!separate_reads(argv[2], argv[3], barcode_len, barcode_stat) )
    {
        return false;
    }
    cout << "   Total time consumed for selecting cells with >= " << effective_cutoff << " read pairs: " 
         << (float)(clock()-tbeg)/CLOCKS_PER_SEC << " seconds.\n" << endl;
    //
    cout << "\n  Total time consumed overall: " 
         << (float)(clock()-tbeg)/CLOCKS_PER_SEC << " seconds.\n" << endl;
    
    return 0;
}
//
bool separate_reads(char* fastqfilename1, char* fastqfilename2, int barcode_len, map<string, unsigned long> barcode_stat)
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
    //
    cout << "\n  Separating reads... " << endl;        
    // parse the original fastq.gz file
    unsigned long checki   = 0;
    unsigned long goocki   = 0;  
    unsigned long goordc   = 0;    
    unsigned long goordc_t = 0;
    int count              = 0;
    int ipart              = 0;
    map<string, string>        barcodebox;
    map<string, unsigned long> barcodereadcnt;    
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
            if(checki%1000000 == 0 && checki!=0)            
            cout << "   Info: " << checki << " read pairs checked.." << endl;
            if(barcodereadcnt.size()%10000 == 0 && barcodereadcnt.size()!=0) 
            cout << "   Info: " << barcodereadcnt.size() << " barcodes now.."  << endl;
            string seq1(""); getline(fp1, seq1);
            string plu1(""); getline(fp1, plu1); 
            string qua1(""); getline(fp1, qua1);            
            string this_barcode = seq1.substr(0, barcode_len);
            string tmpfolder = this_barcode;            
            //
            string seq2(""); getline(fp2, seq2);
            string plu2(""); getline(fp2, plu2);
            string qua2(""); getline(fp2, qua2);
            // if the barcode with sufficient reads according to previous stat
            if(barcode_stat.find(this_barcode) != barcode_stat.end())
            {
                map<string, string>::iterator tmpitr = barcodebox.find(this_barcode);            
                if(tmpitr == barcodebox.end() )
                {            
                    // collect R1: note: remove 16 bp plus 6 bp
                    std::stringstream rss1;
                    rss1.str("");
                    rss1 << line1 << " bc:"  << this_barcode << " UMI:" << seq1.substr(16, 6) << endl
                         << seq1.substr(22)  << endl
                         << plu1             << endl
                         << qua1.substr(22)  << endl;
                    // collect R2
                    std::stringstream rss2;
                    rss2.str("");
                    rss2 << line2 << " bc:"  << this_barcode << endl
                         << seq2  << endl
                         << plu2  << endl
                         << qua2  << endl;                
                    string tmpfiles = "";
                    barcodebox.insert(std::pair<string, string>(this_barcode, tmpfiles));
                    barcodereadcnt.insert(std::pair<string, unsigned long>(this_barcode, 1));                
                    //
                    READ tmpREAD;
                    tmpREAD.R1 = rss1.str();
                    tmpREAD.R2 = rss2.str();
                    vector<READ> tmpVec;
                    tmpVec.push_back(tmpREAD);
                    myBCReads.insert(std::pair<string, vector<READ> >(this_barcode, tmpVec));
                    goocki ++;
                    goordc ++;
                }else
                {
                    string tmpfiles = (*tmpitr).second;
                    map<string, unsigned long>::iterator bcRcntitr = barcodereadcnt.find(this_barcode);
                    assert( bcRcntitr != barcodereadcnt.end() );
                    (*bcRcntitr).second += 1;
                    vector<string> filesinfo = split_string(tmpfiles, '#');
                    // collect R1
                    std::stringstream rss1;
                    rss1.str("");              
                    rss1 << line1 << " bc:"  << this_barcode << " UMI:" << seq1.substr(16, 6) << endl
                         << seq1.substr(22)  << endl
                         << plu1             << endl
                         << qua1.substr(22)  << endl;
                    // collect R2
                    std::stringstream rss2;
                    rss2.str("");                  
                    rss2 << line2 << " bc:"  << this_barcode << endl
                         << seq2  << endl
                         << plu2  << endl
                         << qua2  << endl;
                    //
                    READ tmpREAD;
                    tmpREAD.R1 = rss1.str();
                    tmpREAD.R2 = rss2.str();
                    myBCitr    = myBCReads.find(this_barcode);
                    (*myBCitr).second.push_back(tmpREAD);
                    goocki ++;        
                    goordc ++;                  
                }
            }
        }
        // last output
        if(goocki%min_readpair == 0 && goocki > 0)
        {                      
            ipart ++;
            cout << "   Info: output reads for part " << ipart << " with " << myBCReads.size() << " barcodes: "
                 << goordc << " read pairs. "       << endl;
            goocki    = 0;
            goordc_t += goordc;                      
            goordc    = 0;                 
            if(!output_tmp_reads(&myBCReads, ipart))
            {
                cout << "   Error: output failed. " << endl;
            } 
            myBCReads.clear(); 
            barcodebox.clear();
            barcodereadcnt.clear();  
        }
    }
    // last output
    if(myBCReads.size() > 0)
    {
        ipart ++;                             
        cout << "   Info: output reads for part " << ipart << " with " << myBCReads.size() << " barcodes: "
             << goordc << " read pairs. "       << endl;
        goordc_t += goordc;              
        if(!output_tmp_reads(&myBCReads, ipart))
        {
            cout << "   Error: output failed. " << endl;
        }
    }
    fp1.close();
    fp2.close();
    // close files
    cout << "   Info: "  << checki << " read pairs checked: " << goordc_t << " selected. " << endl;
    //
    // cout << "  Time consumed in sampling : " << (float)(clock()-tbeg)/CLOCKS_PER_SEC << " seconds." << endl;
    return true;
}
//
bool output_tmp_reads(map<string, vector<READ> >* myBCReads, int ipart)
{
    //with current code, reaching here means barcode are good; so no need checking as it is partial.    
    // ipart = checki % 1000000
    unsigned long effectiveBC = 0;   // good barcode with enough read pairs
    unsigned long effectiveRD = 0;   // total number of read pairs corresponding to good barcodes
    map<string, vector<READ> >::iterator bcitr; // barcode separated reads

    map<string, vector<READ> >::iterator bcitr_end;
    bcitr     = (*myBCReads).begin();
    bcitr_end = (*myBCReads).end();
    while(bcitr != bcitr_end)
    {
        string this_barcode       = (*bcitr).first;  // current barcode
        vector<READ> this_readset = (*bcitr).second; // read set for current barcode
        effectiveBC ++;
        effectiveRD += this_readset.size();
        // create new folder if does not exist
        if(!create_folder(output_folder+"/"+this_barcode))
        {
            return false;
        }        
        // open files under folder
        stringstream sepfile1;
        sepfile1.str("");
        sepfile1 << output_folder+"/"+ this_barcode << "/part" << ipart << "_" << this_barcode << "_R1.fq.gz";
        ogzstream* of1 = new ogzstream((sepfile1.str()).c_str(), ios::out);   
        if(!(*of1).good())
        {
             cout << "Cannot open new file " << sepfile1.str() << " to write R1 reads. Exited.\n";
             return false;
        }
        //
        stringstream sepfile2;
        sepfile2.str("");
        sepfile2 << output_folder+"/"+ this_barcode << "/part" << ipart << "_" << this_barcode << "_R2.fq.gz";
        ogzstream* of2 = new ogzstream((sepfile2.str()).c_str(), ios::out);   
        if(!(*of2).good())
        {
             cout << "Cannot open new file " << sepfile2.str() << " to write R2 reads. Exited.\n";
             return false;
        }
        // output reads to partI_*.gz
        vector<READ>::iterator vitr;
        vector<READ>::iterator vitr_end;
        vitr     = this_readset.begin();
        vitr_end = this_readset.end();
        while(vitr != vitr_end)
        {
            READ tmpREAD = *vitr;
            (*of1) << tmpREAD.R1;
            (*of2) << tmpREAD.R2;
            vitr ++;
        }
        //
        (*of1).close();
        (*of2).close();
        // delete reads that have been output to files
        (*myBCReads).erase(bcitr ++);
    }
    //
    return true;
}
//
bool get_barcode_stat(char* fastqfilename1, int barcode_len, int effective_cutoff, map<string, unsigned long>* barcode_stat)
{
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
    map<string, unsigned long> barcode_stat_tmp;
    while(fp1.good())
    {
        string line1("");
        getline(fp1, line1);
        if(line1.size()==0) continue;
        if(line1.substr(0, 1).compare("@")==0) // first line with @XXX is the read id
        {
            // output current read
            string seq1(""); getline(fp1, seq1);
            string plu1(""); getline(fp1, plu1); 
            string qua1(""); getline(fp1, qua1);            
            string this_barcode = seq1.substr(0, barcode_len);
            
            if(barcode_stat_tmp.find(this_barcode) == barcode_stat_tmp.end())
            {
                barcode_stat_tmp.insert(std::pair<string, unsigned long>(this_barcode, 1));
            }else
            {
                barcode_stat_tmp[this_barcode] += 1;
            }
        }
    }
    fp1.close();
    // output
    string ofilename = "asCellseparator_intermediate_raw_barcode_stat.txt";
    ofstream ofp;
    ofp.open(ofilename.c_str(), ios::out);
    if(!ofp.good())
    {
        cout << "   Error: cannot open file to write barcode stat: " << ofilename << endl;
        return false;
    }
    unsigned long total_barco   = barcode_stat_tmp.size();    
    unsigned long total_reads   = 0;
    unsigned long total_barco_o = 0;    
    unsigned long total_reads_o = 0;    
    map<string, unsigned long>::iterator bcsitr;
    map<string, unsigned long>::iterator bcsitr_end;
    bcsitr     = barcode_stat_tmp.begin();
    bcsitr_end = barcode_stat_tmp.end();
    while(bcsitr != bcsitr_end)
    {
        ofp << (*bcsitr).first << "\t" << (*bcsitr).second << endl;
        total_reads += (*bcsitr).second;
        if((*bcsitr).second >= effective_cutoff)
        {
            total_barco_o ++;
            total_reads_o += (*bcsitr).second;
            (*barcode_stat).insert(std::pair<string, unsigned long>((*bcsitr).first, (*bcsitr).second));
        }
        bcsitr ++;
    }
    cout << "   Info: total number of barcodes observed in file1: " << total_barco              << endl;    
    cout << "   Info: total number of reads    observed in file1: " << total_reads              << endl;
    cout << "         among the above, good with read pairs >= "    << effective_cutoff << ": " << endl;
    cout << "               number of barcodes observed in file1: " << total_barco_o            << endl;    
    cout << "               number of reads    observed in file1: " << total_reads_o            << endl;   
    ofp.close();
    //
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

