/*  Funciton: find average sequencing depth in genome-wide sliding windows.
    #
    Hequan Sun, MPIPZ
    Email:sunhequan@gmail.com/sun@mpipz.mpg.de    
*/
#include             <iostream>
#include              <fstream>
#include              <sstream>
#include               <vector>
#include               <string>
#include                  <map>
#include             <assert.h>
#include               <math.h>
// for reading gzipped files
#include   "./gzlib/gzstream.h"
#include               <zlib.h>
// file, dir
#include             <dirent.h>
#include           <sys/stat.h>
#include          <sys/types.h>
#include             <unistd.h>
#include       "split_string.h"
using namespace std;
//
double calculate_windepth(map<unsigned long, double> posDepth, double avgDepth, double* avgDepthReal);
bool update_windepth(map<unsigned long, double>* posDepth, unsigned long winsta);
bool merge_windows(string cnv_raw_file);
//
int main(int argc, char* argv[])
{
    if(argc != 7 && argc != 2)
    {
        cout << "\n  Pls try either function 1: calculate " << endl;
        cout << "    CNV_HQ_v3 1.chrsizes.txt 2.depth_table.txt 3.winSize 4.stepSize 5.basisColumn 6.sampleAvg\n" 
             << endl;
        cout << "       1. gives contigs/chr sizes. "       << endl;
        cout << "       2. gives depth table of multiple samples, each line: chr pos depth1 depth2 ..."        << endl;
        cout << "       3. gives window size "              << endl;
        cout << "       4. gives step   size "              << endl;
      //cout << "       5. gives maximum col-ids as basis for calculating depth-fold change in other samples"  << endl;
        cout << "       5. always be 1."                    << endl;
        cout << "       6. lists the avg depths of all samples (from config.log given by shore consensus)"     << endl;
        cout << endl;
        cout << "    Example: CNV_HQ depth_table.txt 500 100 2 80,60,80,70,80,90,90,75 "                       << endl
             << "              this predicts depth-fold-changes of in genome-wide sliding windows "            << endl
             << "                within windows of 500 bp sliding at 100 bp \n"                                << endl;
        cout << "    Or function 2: merge (need your selection on calculate-output), "                         << endl;
        cout << "\n    CNV_HQ_v3 cnv_hq_mut_wt_diff30percent_upfc.txt\n"                                       << endl;
        //
        return 1;
    }
    //  function 1: calculate window depths
    if(argc == 7)
    {
        // step 0. get inputs
        cout << "   Info: get option settings..." << endl;
        string chrsizesfile        = (string)argv[1];
        string depthfile           = (string)argv[2]; // this is a table with sample-depths in columns
        unsigned long winSize      = strtoul(argv[3], NULL, 0);
        unsigned long step         = strtoul(argv[4], NULL, 0);
        int basis                  = atoi(argv[5]);   // 
        basis                      = 1; // 20200507: always use the first "basis" columns for calculating depth-basis
        vector<string> givendepstr = split_string((string)argv[6], ',');
        vector<double> givendep;
        vector<string>::iterator ditr;
        vector<string>::iterator ditr_end;
        ditr     = givendepstr.begin();
        ditr_end = givendepstr.end();
        while(ditr != ditr_end)
        {
            string tmp = *ditr;
            givendep.push_back(atof(tmp.c_str()));
            ditr ++;
        }
        cout << "   Info: get option settings done." << endl;
        //    
        // step 1: read chr sizes
        cout << "   Info: read chr sizes..." << endl;
        map<string, unsigned long> chrSizes;
        ifstream csfp;
        csfp.open(chrsizesfile.c_str(), ios::in);
        if(!csfp.good())
        {
            cout << "   Error: cannot open file " << chrsizesfile << endl;
            exit(1);
        }
        while(csfp.good())
        {
           string line("");
           getline(csfp, line);
           if(line.size() == 0) continue;
           if(line[0]   == '#') continue;
           vector<string> lineinfo = split_string(line, '\t');  
           if(lineinfo.size()<2) continue;
           string chr = lineinfo[0];
           unsigned long size = strtoul(lineinfo[1].c_str(), NULL, 0);
           if(chrSizes.find(chr) == chrSizes.end())
           {
               chrSizes.insert(std::pair<string, unsigned long>(chr, size));
           }
           else
           {
               cout << "   Warning: repeated chr id found." << endl;
           }
        }
        csfp.close();
        cout << "   Info: read chr sizes done." << endl;
        // prepare output file
        std::stringstream cnvoutfile;
        cnvoutfile.str("");
        cnvoutfile << "cnv_winsize" << winSize << "_step" << step << "_hq.txt";
        ofstream cnvofp;
        cnvofp.open(cnvoutfile.str().c_str(), ios::out);
        // # header: ctg	wstart	wend	ignore	normalized_avg_dep	raw_avg_dep	ctg_size
        // step 2: read and calculate window depth info
        // initialize a map for collecting sample-wise depth in a tmp window
        map<int, map<unsigned long, double> > tmpdepth; // map<sample-id, map<pos, depth> >
        for (int i = 0; i < givendepstr.size(); i ++)
        {
            map<unsigned long, double> tmp;
            tmpdepth.insert(std::pair<int, map<unsigned long, double> >(i, tmp));
        }
        ifstream difp;
        difp.open(depthfile.c_str(), ios::in);
        unsigned long winsta = 1;
        unsigned long winend = winsta + winSize - 1;
        string last_chr("");
        string this_chr("");
        bool firstline  = true; 
        int sample_size = 0;
        double baseDep  = 0;
        cout << "   Info: calculating average seq depth in sliding windows..." << endl;
        while(difp.good())
        {
            string line("");
            getline(difp, line);
            if(line.size() == 0) continue;
            if(line[0]   == '#') continue;            
            // flattened_line_0    4407    86  86  95  83  71  75  65  72
            vector<string> lineinfo = split_string(line, '\t');
            this_chr          = lineinfo[0];
            if(chrSizes.find(this_chr) == chrSizes.end())
            {
                cout << "   Warning: unexpected chr ids in your depth table. Skiped: " << this_chr << endl;
            }
            sample_size = lineinfo.size()-2;
            unsigned long pos = strtoul(lineinfo[1].c_str(), NULL, 0);
            if(firstline) 
            {
                firstline = false;
                last_chr  = lineinfo[0];
            }
            // meet new chr or meet pos out of window for the same chr
            if( (this_chr.compare(last_chr) != 0) ||
                (pos > winend) ) 
            {
                if(winend > chrSizes[last_chr])
                {
                    winend = chrSizes[last_chr];
                }
                // output current info
                // get base average with the first "basis" samples
                baseDep = 0.000000001;
                double avgDepthReal;
                for (int ii = 0; ii < basis; ii ++)
                {
                    double tmpsam = calculate_windepth(tmpdepth[ii], givendep[ii], &avgDepthReal);
                    baseDep += tmpsam;
                }   
                baseDep = baseDep / (basis+0.0000001);
                if(baseDep < 1) baseDep += 1; // caution 0
                // get fold-change in other samples
                // if(baseDep>=30) // caution minimum coverage applied; turned off on 2019-11-28
                std::stringstream sscov;
                sscov.str("");               
                for (int ii = 0; ii < sample_size; ii ++)
                {
                    double tmpsam = calculate_windepth(tmpdepth[ii], givendep[ii], &avgDepthReal);
                    if(tmpsam == 0) tmpsam += 0.01*baseDep;
                    // common
                    if(ii == 0)
                    {
                        cnvofp << last_chr << "\t" << winsta << "\t" << winend << "\t";
                    }
                    // sample
                    cnvofp << roundf(tmpsam / baseDep * 100) / 100 << "\t"; 
                    sscov  << roundf(tmpsam / 1.00000 * 1)   / 1   << "\t" << avgDepthReal << "\t"; 
                    // real depth scaled against sample-average, and non-scaled
                }
                cnvofp << sscov.str() << chrSizes[last_chr] << endl; // real depth scaled against sample-average
                // update next chr
                if(this_chr.compare(last_chr) != 0 )
                {
                    last_chr = this_chr;
                    winsta = 1;
                    winend = winsta + winSize - 1;
                    // remove depth info not in new chr
                    for (int ii = 0; ii < sample_size; ii ++)
                    {
                        tmpdepth[ii].clear();
                    }
                }
                else
                if(pos > winend)
                {
                    winsta = winsta + step; 
                    winend = winsta + winSize - 1; 
                    if(winend > chrSizes[this_chr])
                    {
                        winend = chrSizes[this_chr];
                    }
                    // remove depth info not in next window but same chr
                    for (int ii = 0; ii < sample_size; ii ++)
                    {
                        update_windepth(&tmpdepth[ii], winsta);
                    }
                }
            }
            // collect sample-wise position depth: 2-10 means samples CDAB/MNOP in my apple case
            for (int jj = 2; jj < lineinfo.size(); jj ++)
            {
                double idepth = atof(lineinfo[jj].c_str());
                tmpdepth[jj - 2].insert(std::pair<unsigned long, double>(pos, idepth));
            }      
        }
        // last window
        cout << "   Info: checking last window..." << endl;
        if(baseDep == 0)
        {
            baseDep = 0.000000001;
            double avgDepthReal;
            for (int ii = 0; ii < basis; ii ++)
            {
                double tmpsam = calculate_windepth(tmpdepth[ii], givendep[ii], &avgDepthReal);
                baseDep += tmpsam;
            }
            baseDep = baseDep / (basis+0.0000001);
            if(baseDep < 1) baseDep += 1; // caution 0
        }
        if(winend > chrSizes[last_chr])
        {
            winend = chrSizes[last_chr];
            cout << "   Info: winend changed to " << chrSizes[last_chr] << endl;
        }
        std::stringstream sscov;
        sscov.str("");          
        for (int ii = 0; ii < sample_size; ii ++)
        {
            double avgDepthReal;
            double tmpsam = calculate_windepth(tmpdepth[ii], givendep[ii], &avgDepthReal);
            if(tmpsam == 0) tmpsam += 0.01*baseDep;            
            // common
            if(ii == 0)
            {
                cnvofp << last_chr << "\t" << winsta << "\t" << winend << "\t";
            }
            // sample
            cnvofp << roundf(tmpsam / baseDep * 100) / 100 << "\t"; 
            sscov  << roundf(tmpsam / 1.00000 * 1  ) / 1   << "\t" << avgDepthReal << "\t"; 
            // real depth scaled against sample-average, and non-scaled           
        }
        cnvofp << sscov.str() << chrSizes[last_chr] << endl;  // real depth scaled against sample-average
        cout << "   Info: checking last window done." << endl;
        cout << "   Info: calculating average seq depth in sliding windows done." << endl;
        //
        difp.close();
        cnvofp.close();
        //
        cout << "   Info: average depth in sliding windows have been collected in " << cnvoutfile.str() << endl;
    }
    //  function 2: merge overlapping windows on the same contigs
    if(argc == 2)
    {
        if(!merge_windows((string)(argv[1])))
        {
            return 1;
        }
    } 
    return 0;
}
// sample-wise
double calculate_windepth(map<unsigned long, double> posDepth, double avgDepth, double* avgDepthReal)
{
    // input:
    // posDepth : individual depths at positions within the window
    // avgDepth : sample-wise genome-wide average depth from input
    // output
    // iwinDepth: averge depth within the given window, and normalized to the given genome-wide average depth     
    double iwinDepth = 0;
    //
    map<unsigned long, double>::iterator pitr;
    map<unsigned long, double>::iterator pitr_end;
    pitr     = posDepth.begin();
    pitr_end = posDepth.end();
    while(pitr != pitr_end)
    {
        iwinDepth += (*pitr).second;
        pitr ++;
    } 
    // real value
    *avgDepthReal = iwinDepth / posDepth.size();
    //
    iwinDepth = iwinDepth / posDepth.size() / avgDepth * 100; // "normalized" with avgDepth * 100
    //
    return iwinDepth;
}
//
bool update_windepth(map<unsigned long, double>* posDepth, unsigned long winsta)
{
    // remove info before winsta
    map<unsigned long, double>::iterator pitr;
    map<unsigned long, double>::iterator pitr_end;
    pitr     = (*posDepth).begin();
    pitr_end = (*posDepth).end();
    while(pitr != pitr_end)
    {
        if((*pitr).first < winsta)
        {
            (*posDepth).erase(pitr++);
        }
        else pitr ++;
    } 
    return true; 
}
//
bool merge_windows(string cnv_raw_file)
{
    // prepare output
    ofstream upfp;
    string ofilename = cnv_raw_file + ".merged";
    upfp.open(ofilename.c_str(), ios::out); // merging selected sliding windows
    if(!upfp.good())
    {
        cout << "   Error: cannot open file " << "cnv_hq_window_merged_upMUT_final.txt" << endl;
        return false;
    }
    // raw fold change file from previous step
    ifstream ifp;
    ifp.open(cnv_raw_file.c_str(), ios::in);
    if(!ifp.good())
    {
        cout << "   Error: cannot find file " << cnv_raw_file << endl;
        
    }   
    string this_chr("");
    string last_chr("");
    bool firstline = true; 
    vector<double> linedepthsum;
    int sample_size = 0;
    int line_merged = 0;
    unsigned long regionsta;
    unsigned long regionend;
    while(ifp.good())
    {
        string line("");
        getline(ifp, line);
        if(line.size()==0 || line[0]=='#') continue;
        // flattened_line_0    1751-1950   1.12    0.88    0.68    0.77    1.08    0.95    0.91    0.75
        vector<string> lineinfo = split_string(line, '\t');
        this_chr                = lineinfo[0];                       
        vector<string> posinfo  = split_string(lineinfo[1], '-'); 
        unsigned long this_sta  = strtoul(posinfo[0].c_str(), NULL, 0);
        unsigned long this_end  = strtoul(posinfo[1].c_str(), NULL, 0);
        
        if(firstline)
        {
            firstline   = false;
            last_chr    = this_chr;
            sample_size = lineinfo.size()-2; // = 8
            for (int ii = 0; ii < sample_size; ii ++)
            {
                double tmpd = atof(lineinfo[ii+2].c_str()); // col: 2-10
                linedepthsum.push_back(tmpd);
            }
            line_merged = 1;
            regionsta = this_sta;
            regionend = this_end;
        }
        else
        {
            if(this_chr.compare(last_chr) != 0) // meet new chr
            {
                // output
                upfp << last_chr << "\t" << regionsta << "-" << regionend << "\t";
                for(int ii = 0; ii < sample_size; ii ++)
                {
                    if(ii < sample_size - 1)
                    {
                        upfp << roundf( linedepthsum[ii] / line_merged *100 ) / 100 << "\t";
                    }
                    else
                    {
                        upfp << roundf( linedepthsum[ii] / line_merged *100 ) / 100 << "\tavg_of_" 
                             << line_merged << "_lines" << endl;
                    }
                }
                // new chr
                last_chr    = this_chr;  
                for (int ii = 0; ii < sample_size; ii ++)
                {
                    double tmpd = atof(lineinfo[ii+2].c_str()); // col: 2-10
                    linedepthsum[ii] = tmpd;
                }
                line_merged = 1;
                regionsta = this_sta;
                regionend = this_end;                              
            }else // meet new region of the same chr
            if(this_sta > regionend)
            {
                // output
                upfp << last_chr << "\t" << regionsta << "-" << regionend << "\t";
                for(int ii = 0; ii < sample_size; ii ++)
                {
                    if(ii < sample_size - 1)
                    {
                        upfp << roundf( linedepthsum[ii] / line_merged *100 ) / 100 << "\t";
                    }
                    else
                    {
                        upfp << roundf( linedepthsum[ii] / line_merged *100 ) / 100 << "\tavg_of_" 
                             << line_merged << "_lines" << endl;
                    }
                }
                // new chr
                last_chr    = this_chr;  
                for (int ii = 0; ii < sample_size; ii ++)
                {
                    double tmpd = atof(lineinfo[ii+2].c_str()); // col: 2-10
                    linedepthsum[ii] = tmpd;
                }
                line_merged = 1;
                regionsta = this_sta;
                regionend = this_end;                   
            }else
            {
                for (int ii = 0; ii < sample_size; ii ++)
                {
                    double tmpd = atof(lineinfo[ii+2].c_str()); // col: 2-10
                    linedepthsum[ii] += tmpd; // take the sum of the same region, later make average
                }
                line_merged ++;
                regionend = this_end; 
                last_chr    = this_chr;                                 
            }
        }
    }  
    // output last region
    // output
    upfp << last_chr << "\t" << regionsta << "-" << regionend << "\t";
    for(int ii = 0; ii < sample_size; ii ++)
    {
        if(ii < sample_size - 1)
        {
            upfp << roundf( linedepthsum[ii] / line_merged *100 ) / 100 << "\t";
        }
        else
        {
            upfp << roundf( linedepthsum[ii] / line_merged *100 ) / 100 << "\tavg_of_" 
                 << line_merged << "_lines" << endl;
        }
    }          
    //
    ifp.close();
    upfp.close();
    return true;
}



