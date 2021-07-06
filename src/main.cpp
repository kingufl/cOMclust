#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>
#include <set>
#include <ctime>
#include <vector>
#include <stdlib.h>
#include <typeinfo>
#include <string>
#include "KDTree.hpp"

#define NDEBUG
#include <cassert>

#define Pi 3.14159265
#define pi 3.14156

using namespace std;

inline double dist2(const vector<float> &a, vector<float> &b) {
    double distc = 0;
    for (size_t i = 0; i < a.size(); i++) {
        double di = a.at(i) - b.at(i);
        distc += di * di;
    }
    return sqrt(distc);
}

vector<vector<float> > data;

struct trip{

    int first;
    int second;
    int third;

};


pair<int,int> count21=make_pair(0,0), count31=make_pair(0,0), count32=make_pair(0,0);

vector < vector <float> > corrected_rmaps;
vector<string> rmap_belong;


vector < int > rmaps_to_correct;

vector< vector<float> > rmap,rmap_rev;

vector< vector< int> > omrels;

vector< vector<int> > rmap_quantized;

vector < vector <int> > rmap_store;

vector < set <int> > rmap_relations;

vector < vector <int> > rmap_relations_filtered;

vector < vector <int> > rmap_relations_final;

vector <vector < vector < pair <int, int> > > > rmap_alignments;

vector <vector < vector < trip > > > multi_align_grid;

vector< string > rmap_name;

vector < vector < pair <int, int> > > kmer_store;

vector < vector < string > > concensus_res;

int buket_size=2;
map<string,int> nameid;


map < string, int > kmer_map;

map < int, string > reverse_map;

vector < pair <int, int> > alignment;

int max_alignments=30; // Number of related Rmaps used for error correcting an Rmap. Use 50 for improved error correction with an increase in run time.

int k_size=4;

int m_size=3;

int min_frags=15;

int align_size=3;

float threshold=0.85;

int min_alignments = 3;

int min_concensus = 2; // min consensus needed on a fragment to error correct it. Increase to 3 or 4 for greater confidence. Should increase max_alignments to 50 as well.

float bias_for_single=0.1;

int start_r,end_r,diff;

string nu;

int rmap_count;

float find_score( float a, float b){

    if (b>a)
        return(a/b);

    else
        return(b/a);
}

struct triplet{

    int pos1;
    int rmap_2;
    int pos2;

};


vector< vector < triplet > > complete_relations;

vector < pair<int ,int> > optimized_overlap_alignment(vector< float >& rmap1, vector< float >& rmap2, int p1, int p2);


float unit_bias=0.2;

trip make_trip(int a, int b, int c){

    trip t;
    t.first=a;
    t.second=b;
    t.third=c;
    return(t);
}

int check_relation(int rmap1, int pos1, int rmap2, int pos2){

//    return 1;

    int flag=1;

    //cout << endl << endl;

    //cout << "Aligning " << rmap1 << " with " << rmap2 <<endl ;
    float total_sc=0;

    while((pos1<rmap[rmap1].size()-1) && (pos2<rmap[rmap2].size()-1)){

        float sc=0;

        float score=0, new_score=0;

        int cp1=pos1, cp2=pos2;

        score= find_score(rmap[rmap1][pos1],rmap[rmap2][pos2])+bias_for_single; // 1-1
        sc=score/2;

//        if(cp1+1<rmap[rmap1].size() && cp2+1<rmap[rmap2].size())
//            score= score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][cp1+1],rmap[rmap2][cp2+1]);


        //cout << "compare " << rmap[rmap1][pos1] << "  AND  " << rmap[rmap2][pos2] << " Score : " << score << endl;
        if(pos1+1<rmap[rmap1].size()){

            new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1],rmap[rmap2][pos2]); //2-1

            if(pos1+2<rmap[rmap1].size() && pos2+1<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+2],rmap[rmap2][pos2+1]);

            //cout << "compare " << rmap[rmap1][pos1]<< "+"<< rmap[rmap1][pos1+1] << "  AND  " << rmap[rmap2][pos2] << " Score : " << new_score<< endl;

            if (new_score>score){
                    score=new_score;
                    cp1=pos1+1;
                    cp2=pos2;

                    sc=score/3;


            }
        }


        if(pos1+2<rmap[rmap1].size())
        new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1]+rmap[rmap1][pos1+2],rmap[rmap2][pos2]); //3-1

        if(pos1+3<rmap[rmap1].size() && pos2+1<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+3],rmap[rmap2][pos2+1]);

        //cout << "compare " << rmap[rmap1][pos1]<< "+"<< rmap[rmap1][pos1+1] << "+"<< rmap[rmap1][pos1+2] << "  AND  " << rmap[rmap2][pos2] << " Score : " << new_score<< endl;


        if (new_score>score){
                score=new_score;
                cp1=pos1+2;
                cp2=pos2;
                sc=score/4;
        }

        if(pos2+1<rmap[rmap2].size())
        new_score= find_score(rmap[rmap1][pos1],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]); //1-2

        if(pos1+1<rmap[rmap1].size() && pos2+2<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+1],rmap[rmap2][pos2+2]);

        if (new_score>score){
                score=new_score;
                cp1=pos1;
                cp2=pos2+1;
                sc=score/3;
        }

        if(pos2+2<rmap[rmap2].size())
        new_score= find_score(rmap[rmap1][pos1],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]+rmap[rmap2][pos2+2]); //1-3

        if(pos1+1<rmap[rmap1].size() && pos2+3<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+1],rmap[rmap2][pos2+3]);

        if (new_score>score){
                score=new_score;
                cp1=pos1;
                cp2=pos2+2;
                sc=score/4;
        }

        if((pos1+1<rmap[rmap1].size()) && (pos2+1<rmap[rmap2].size()))
        new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]);//

        if(pos1+2<rmap[rmap1].size() && pos2+2<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+2],rmap[rmap2][pos2+2]);

        if (new_score>score){
                    score=new_score;
                    cp1=pos1+1;
                    cp2=pos2+1;
                    sc=score/4;
        }

        if((pos1+1<rmap[rmap1].size()) && (pos2+2<rmap[rmap2].size()))
            new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]+rmap[rmap2][pos2+2]);

        if(pos1+2<rmap[rmap1].size() && pos2+3<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+2],rmap[rmap2][pos2+3]);

        if (new_score>score){
                    score=new_score;
                    cp1=pos1+1;
                    cp2=pos2+2;
                    sc=score/5;
        }

        if((pos1+2<rmap[rmap1].size()) && (pos2+1<rmap[rmap2].size()))
            new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1]+rmap[rmap1][pos1+2],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]);

        if(pos1+3<rmap[rmap1].size() && pos2+2<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+3],rmap[rmap2][pos2+2]);

        if (new_score>score){
                    score=new_score;
                    cp1=pos1+2;
                    cp2=pos2+1;
                    sc=score/5;
        }


        if((pos1+2<rmap[rmap1].size()) && (pos2+2<rmap[rmap2].size()))
            new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1]+rmap[rmap1][pos1+2],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]+rmap[rmap2][pos2+2]);

        if(pos1+3<rmap[rmap1].size() && pos2+3<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+3],rmap[rmap2][pos2+3]);

        if (new_score>score){
                    score=new_score;
                    cp1=pos1+2;
                    cp2=pos2+2;
                    sc=score/6;
        }

        if((pos1+3<rmap[rmap1].size()) && (pos2+2<rmap[rmap2].size()))
            new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1]+rmap[rmap1][pos1+2]+rmap[rmap1][pos1+3],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]+rmap[rmap2][pos2+2]);

        if(pos1+4<rmap[rmap1].size() && pos2+3<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+4],rmap[rmap2][pos2+3]);

        if (new_score>score){
                    score=new_score;
                    cp1=pos1+3;
                    cp2=pos2+2;
                    sc=score/7;
        }

        if((pos1+2<rmap[rmap1].size()) && (pos2+3<rmap[rmap2].size()))
            new_score= find_score(rmap[rmap1][pos1]+rmap[rmap1][pos1+1]+rmap[rmap1][pos1+2],rmap[rmap2][pos2]+rmap[rmap2][pos2+1]+rmap[rmap2][pos2+2]+rmap[rmap2][pos2+3]);

        if(pos1+3<rmap[rmap1].size() && pos2+4<rmap[rmap2].size())
                new_score= new_score*(1-unit_bias)+ unit_bias*find_score(rmap[rmap1][pos1+3],rmap[rmap2][pos2+4]);

        if (new_score>score){
                    score=new_score;
                    cp1=pos1+2;
                    cp2=pos2+3;
                    sc=score/7;
        }



        if (score>threshold){
            total_sc+=sc;

//            for(int i=pos1; i<=cp1; i++)
//                cout << i << " (" << rmap[rmap1][i] << ")   ";
//
//            cout << " aligns with ";
//            for(int i=pos2; i<=cp2; i++)
//                cout << i << " (" << rmap[rmap2][i] << ")   ";
//
//            cout << endl;
        }

        else{

//            cout << "Maps do not align ";
            flag=0;

            break;

        }

        pos1=cp1+1;
        pos2=cp2+1;

    }
    if(total_sc<2)
        flag=0;
    //if((pos1>=rmap[rmap1].size()) || (pos2>=rmap[rmap2].size()))
        //cout <<endl << " Alignment successful";
    //cout << "  FLAG: " << flag;

    if(flag==1){
        //cout<< endl << "Alignment success" << endl;
    }
        return(flag);
}



void read_data( char* filename){

    rmap_count=0;

    cout<< "Reading rmaps....."<<endl;

    vector<string> rnames,rrnames;

    ifstream file(filename, ifstream::in);
	if (!file.is_open()) {
		cout << "Error opening file" << endl;
		exit(1);
	}

	int total_rmaps=0;
//	cout << "Number of Rmaps ? ";

	//cin >> total_rmaps;

	vector<vector<float> > data_r;
	vector<string> rmap_belong_r;




    int kmercnt=0;
	string str,str2;
    int lastname=0;
	while(getline(file, str))
    {
        istringstream s(str);

        vector < float > v1;
        vector < int > v2;

        s >> str2;

        string revst = str2+"_r";
        total_rmaps++;

        getline(file, str);

        //cout << rmap_count << " ";

//        float first_num;
//        ss >> first_num;
        istringstream ss(str);

        string dummy;
        ss>>dummy>>dummy;

        float num;
        while(ss >> num)
        {
        //cout << num << "    ";
         v1.push_back(num);
         v2.push_back(round(num/buket_size));
        }

        if(v1.size()<4)
            continue;

        rnames.push_back(str2);

        vector<float> newkmer(4,0);
        for(int ii=0;ii<v1.size()-4+1;ii++){
            int index=0;
            for(int jj=ii;jj<4+ii;jj++){
                newkmer[index++]=asinh(v1[jj]);
            }
            data.push_back(newkmer);
            rmap_belong.push_back(str2);
            kmercnt++;
        }

        //cout<<str2<<" "<<"read"<<endl;
        rmap_name.push_back(str2);
        rmap.push_back(v1);

        //Uncomment the following part to increase coverage to 2x. Rmaps are then considered in both forward and reverse direction.
        //This increases memory and running time by 2x, but delivers slightly improved error-correction.

        /*
        reverse(v1.begin(),v1.end());
        rmap_rev.push_back(v1);
        rrnames.push_back(revst);

        for(int ii=0;ii<v1.size()-4+1;ii++){
            int index=0;
            for(int jj=ii;jj<4+ii;jj++){
                newkmer[index++]=asinh(v1[jj]);
            }
            data_r.push_back(newkmer);
            rmap_belong_r.push_back(revst);
            kmercnt++;
        }

        */

        getline(file, str);
        rmap_count++;


    }
    //cout << endl;
    file.close();

    rmap.insert(rmap.end(),rmap_rev.begin(),rmap_rev.end());
    rnames.insert(rnames.end(),rrnames.begin(),rrnames.end());
    data.insert(data.end(),data_r.begin(),data_r.end());
    rmap_belong.insert(rmap_belong.end(),rmap_belong_r.begin(),rmap_belong_r.end());

    for(int i=0;i<rnames.size();i++){
        if(nameid.find(rnames[i])==nameid.end()){
            nameid[rnames[i]]=lastname++;
        }
    }

    cout << "Total rmaps : " <<total_rmaps<<endl<<"Rmaps with minimum fragments : "<< rmap_count<<endl;

    cout<<"total kmers: "<<kmercnt<<endl;

    total_rmaps=rmap_count;

//	rmap_relations.resize(total_rmaps);
//	rmap_relations_filtered.resize(total_rmaps);
	rmap_relations_final.resize(total_rmaps);
	rmap_alignments.resize(total_rmaps);
	multi_align_grid.resize(total_rmaps);
	concensus_res.resize(total_rmaps);

//	complete_relations.resize(total_rmaps);

    cout<<"Rmaps read"<<endl;

}


void build_multi_align_grid(){


    long mag=0;
    for(int i=0; i<rmap_count; i++){


        for(int j=0; j<rmap_alignments[i].size(); j++){

            mag+=rmap_alignments[i][j].size();
            for(int k=1; k<rmap_alignments[i][j].size()-1; k++){

                int ref_in=rmap_alignments[i][j][k].first;
                int tar_in=rmap_alignments[i][j][k].second;

//                if(k==rmap_alignments[i][j].size()-1){
//                    multi_align_grid[i][j][ref_in]=make_pair(ref_gap,tar_gap);
//
//                }

                int ref_gap=rmap_alignments[i][j][k+1].first-rmap_alignments[i][j][k].first;
                int tar_gap=rmap_alignments[i][j][k+1].second-rmap_alignments[i][j][k].second;

                while(ref_in<rmap_alignments[i][j][k+1].first){
                    multi_align_grid[i][j][ref_in]=make_trip(ref_gap,tar_gap,tar_in);
                    ref_in++;
                }


                }

            }

    }

}


void find_concensus(char* ch){



    int lastmapped=0;

    ofstream out;

    string st;
    string str(ch);

    st = "results/concensus"+str+".txt";
    out.open(st.c_str());


    for(int i=0; i<rmap_count; i++){

    rmaps_to_correct.push_back(i);

    if(multi_align_grid[i].size()<min_alignments)
        continue;

    out<<rmap_name[i];

//    cout << multi_align_grid[i][0].size()<<endl;

        for(int k=0; k<multi_align_grid[i][0].size(); k++){


            lastmapped=0;
            map < string, int > mapit;
            map < int, string > rev_mapit;
            vector < int > keep_count(multi_align_grid[i].size(),0);

            for(int j=0; j<multi_align_grid[i].size(); j++){

                if(multi_align_grid[i][j][k].first==0)
                    continue;

                string str;
                ostringstream oss;
                oss<<multi_align_grid[i][j][k].first<<","<<multi_align_grid[i][j][k].second;

                istringstream iss(oss.str());
                iss >> str;

                if(mapit.find(str)==mapit.end()){
                    mapit[str]=lastmapped;
                    rev_mapit[lastmapped]=str;
                    lastmapped++;
                    keep_count[mapit[str]]++;

                }
                else{
                    keep_count[mapit[str]]++;
                }

            }
            int maxel=0;
            for(int x=0; x<keep_count.size(); x++){
                if(keep_count[x]>keep_count[maxel])
                    maxel=x;
            }

            if(keep_count[maxel]<min_concensus){
                concensus_res[i].push_back("0,0");
                out << " (0,0) ";
            }
            else{
                concensus_res[i].push_back(rev_mapit[maxel]);
                out << " (" <<rev_mapit[maxel]<<") ";
            }



        }
        out<<endl;

    }

out.close();

}


bool check(string & str, int x, int y){

                char chara=str[0];
                char charb=str[2];

                int a=chara - '0';
                int b=charb - '0';

                if(a==x && b==y)
                    return true;
                else
                    return false;
}

void correct_rmaps(char* ch){

    ofstream out2;
    ofstream out;

    long tot_frag=0;
    long untouched=0;

    string st;

    string in(ch);

    st = "results/corrections"+in+".txt";

    out2.open(st.c_str());

    st=  "results/corrected_rmaps"+in+".txt";
    out.open(st.c_str());


    corrected_rmaps.resize(rmap_count);

    for(int i=0; i<rmap_count; i++){

//        if(i%1000==0){
//            cout<<i<<endl;
//        }

        if(multi_align_grid[i].size()<min_alignments){
            continue;
            out<<rmap_name[i]<<"_nc"<<endl<<"enzyme enzyme ";

            for(int ii=0;ii<rmap[i].size();ii++){
                out<<rmap[i][ii]<<" ";
            }

            out<<endl<<endl;

            continue;
        }


        int add_err=0;
        int del_err=0;
        int this_rmap=i;
        out<<rmap_name[this_rmap]<<endl<<"enzyme enzyme ";

        bool set_bit2=false;
        bool set_bit3=false;
        bool set_bit4=false;
        float counter;

        tot_frag+=concensus_res[this_rmap].size();

        for(int j=0; j<concensus_res[this_rmap].size(); j++){

            if (set_bit2==true){
                    set_bit2=false;
                    continue;
            }

            else if(set_bit3==true){
                set_bit3=false;
                set_bit2=true;
                continue;
            }

            else if(set_bit4==true){
                set_bit4=false;
                set_bit3=true;
                continue;
            }

            string str = concensus_res[this_rmap][j];

            char chara=str[0];
            char charb=str[2];

            int a=chara - '0';
            int b=charb - '0';

            if(a==1 && b==1){
                float tot11=0;
                counter=0;
                for(int k=0; k<multi_align_grid[this_rmap].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot11+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                    }
                }
                tot11+=rmap[this_rmap][j];
                counter++;

                corrected_rmaps[i].push_back(tot11/counter);
//                out<<"  "<<tot11/counter<<"("<<j<<","<<str<<"correction)";
                out<<"  "<<roundf((tot11/counter) * 1000) / 1000;


            }

            else if(a==1 && b==2){
                del_err++;
                float tot121=0;
                float tot122=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot121+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot122+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                    }

                }
                corrected_rmaps[i].push_back(tot121/counter);
                corrected_rmaps[i].push_back(tot122/counter);
//                out << "  "<<tot121/counter<<"  "<<tot122/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot121/counter) * 1000) / 1000<<"  "<<roundf((tot122/counter) * 1000) / 1000;
            }

            else if(a==2 && b==1){


                count21.first++;
                if(!check(concensus_res[this_rmap][j+1],2,1)){
                    count21.second++;
                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
                    out << "  "<<rmap[this_rmap][j];
                    untouched++;
                    continue;
                }
                add_err++;
                float tot21=0;
                set_bit2=true;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot21+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                    }

                }
                corrected_rmaps[i].push_back(tot21/counter);
//                out << "  "<<tot21/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot21/counter) * 1000) / 1000;

            }

            else if(a==2 && b==3){

            del_err++;

//                count21.first++;
//                if(!check(concensus_res[this_rmap][j+1],2,1)){
//                    count21.second++;
//                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
//                    out << "  "<<rmap[this_rmap][j];
//                    continue;
//                }
                float tot231=0;
                float tot232=0;
                float tot233=0;
                counter=0;
                set_bit2=true;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot231+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot232+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                        tot233+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+2];
                    }

                }
                corrected_rmaps[i].push_back(tot231/counter);
                corrected_rmaps[i].push_back(tot232/counter);
                corrected_rmaps[i].push_back(tot233/counter);
//                out << "  "<<tot21/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot231/counter) * 1000) / 1000<< "  "<<roundf((tot232/counter) * 1000) / 1000<< "  "<<roundf((tot233/counter) * 1000) / 1000;

            }

             else if(a==2 && b==4){
                del_err++;
                del_err++;

//                count21.first++;
//                if(!check(concensus_res[this_rmap][j+1],2,1)){
//                    count21.second++;
//                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
//                    out << "  "<<rmap[this_rmap][j];
//                    continue;
//                }
                float tot241=0;
                float tot242=0;
                float tot243=0;
                float tot244=0;
                set_bit2=true;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot241+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot242+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                        tot243+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+2];
                        tot244+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+3];
                    }

                }
                corrected_rmaps[i].push_back(tot241/counter);
                corrected_rmaps[i].push_back(tot242/counter);
                corrected_rmaps[i].push_back(tot243/counter);
                corrected_rmaps[i].push_back(tot244/counter);
//                out << "  "<<tot21/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot241/counter) * 1000) / 1000<< "  "<<roundf((tot242/counter) * 1000) / 1000<< "  "<<roundf((tot243/counter) * 1000) / 1000<< "  "<<roundf((tot244/counter) * 1000) / 1000;

            }

            else if(a==2 && b==2){

                if(!check(concensus_res[this_rmap][j+1],2,2)){
                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
                    out << "  "<<rmap[this_rmap][j];
                    continue;
                }
                set_bit2=true;
                float tot221=0;
                float tot222=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot221+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot222+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                    }

                }
                corrected_rmaps[i].push_back(tot221/counter);
                corrected_rmaps[i].push_back(tot222/counter);
//                out << "  "<<tot221/counter<<"  "<<tot222/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot221/counter) * 1000) / 1000<<"  "<<roundf((tot222/counter) * 1000) / 1000;

            }

            else if(a==1 && b==3){
                del_err++;
                del_err++;
                float tot131=0;
                float tot132=0;
                float tot133=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot131+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot132+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                        tot133+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+2];
                    }

                }
                corrected_rmaps[i].push_back(tot131/counter);
                corrected_rmaps[i].push_back(tot132/counter);
                corrected_rmaps[i].push_back(tot133/counter);
//                out << "  "<<tot121/counter<<"  "<<tot122/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot131/counter) * 1000) / 1000<<"  "<<roundf((tot132/counter) * 1000) / 1000<< "  " << roundf((tot133/counter) * 1000) / 1000;
            }

            else if(a==1 && b==4){
                del_err++;
                del_err++;
                del_err++;
                float tot141=0;
                float tot142=0;
                float tot143=0;
                float tot144=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot141+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot142+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                        tot143+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+2];
                        tot144+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+3];
                    }

                }
                corrected_rmaps[i].push_back(tot141/counter);
                corrected_rmaps[i].push_back(tot142/counter);
                corrected_rmaps[i].push_back(tot143/counter);
                corrected_rmaps[i].push_back(tot144/counter);
//                out << "  "<<tot121/counter<<"  "<<tot122/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot141/counter) * 1000) / 1000<<"  "<<roundf((tot142/counter) * 1000) / 1000<< "  " << roundf((tot143/counter) * 1000) / 1000<< "  " << roundf((tot144/counter) * 1000) / 1000;
            }

            else if(a==3 && b==3){


                set_bit3=true;
                float tot331=0;
                float tot332=0;
                float tot333=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot331+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot332+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                        tot333+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+2];
                    }

                }
                corrected_rmaps[i].push_back(tot331/counter);
                corrected_rmaps[i].push_back(tot332/counter);
                corrected_rmaps[i].push_back(tot333/counter);
//                out << "  "<<tot331/counter<<"  "<<tot332/counter<<"  "<<tot333/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot331/counter) * 1000) / 1000<<"  "<<roundf((tot332/counter) * 1000) / 1000<<"  "<<roundf((tot333/counter) * 1000) / 1000;

            }

            else if(a==3 && b==1){


                count31.first++;
                if(!(check(concensus_res[this_rmap][j+1],3,1))){
                    count31.second++;
                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
                    out << "  "<<rmap[this_rmap][j];
                    untouched++;
                    continue;
                }
                add_err++;
                add_err++;

                set_bit3=true;
                float tot31=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot31+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                    }

                }
                corrected_rmaps[i].push_back(tot31/counter);
//                out << "  "<<tot321/counter<<"  "<<tot322/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot31/counter) * 1000) / 1000;

            }

            else if(a==3 && b==2){


                count32.first++;
                if(!(check(concensus_res[this_rmap][j+1],3,2))){
                    count32.second++;
                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
                    out << "  "<<rmap[this_rmap][j];
                    untouched++;
                    continue;
                }
                add_err++;
                set_bit3=true;
                float tot321=0;
                float tot322=0;
                counter=0;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot321+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot322+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                    }

                }
                corrected_rmaps[i].push_back(tot321/counter);
                corrected_rmaps[i].push_back(tot322/counter);
//                out << "  "<<tot321/counter<<"  "<<tot322/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot321/counter) * 1000) / 1000<<"  "<<roundf((tot322/counter) * 1000) / 1000;

            }

            else if(a==4 && b==1){

            add_err++;
            add_err++;
            add_err++;
//                count21.first++;
//                if(!check(concensus_res[this_rmap][j+1],2,1)){
//                    count21.second++;
//                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
//                    out << "  "<<rmap[this_rmap][j];
//                    continue;
//                }
                float tot41=0;
                counter=0;
                set_bit4=true;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot41+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];

                    }

                }
                corrected_rmaps[i].push_back(tot41/counter);
//                out << "  "<<tot21/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot41/counter) * 1000) / 1000;

            }

            else if(a==4 && b==3){

            add_err++;
//                count21.first++;
//                if(!check(concensus_res[this_rmap][j+1],2,1)){
//                    count21.second++;
//                    corrected_rmaps[i].push_back(rmap[this_rmap][j]);
//                    out << "  "<<rmap[this_rmap][j];
//                    continue;
//                }
                float tot431=0;
                float tot432=0;
                float tot433=0;
                counter=0;
                set_bit4=true;
                for(int k=0; k<multi_align_grid[rmaps_to_correct[i]].size(); k++){

                    if(multi_align_grid[this_rmap][k][j].first==a && multi_align_grid[this_rmap][k][j].second==b){
                        counter++;
                        tot431+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third];
                        tot432+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+1];
                        tot433+=rmap[rmap_relations_final[this_rmap][k]][multi_align_grid[this_rmap][k][j].third+2];
                    }

                }
                corrected_rmaps[i].push_back(tot431/counter);
                corrected_rmaps[i].push_back(tot432/counter);
                corrected_rmaps[i].push_back(tot433/counter);
//                out << "  "<<tot21/counter<<"("<<j<<","<<str<<"correction)";
                out << "  "<<roundf((tot431/counter) * 1000) / 1000<< "  "<<roundf((tot432/counter) * 1000) / 1000<< "  "<<roundf((tot433/counter) * 100) / 100;

            }



            else if(a==0 && b==0){
                corrected_rmaps[i].push_back(rmap[this_rmap][j]);
                out << "  "<<rmap[this_rmap][j];
                untouched++;
            }

            else{
                corrected_rmaps[i].push_back(rmap[this_rmap][j]);
                out << "  "<<rmap[this_rmap][j];
            }


        }

        out2<< rmap_name[this_rmap]<<"  "<<add_err<<"   "<< del_err<<endl;

        out<<endl<<endl;

    }

    out.close();
    out2.close();

}



int main ( int argc, char *argv[] ) {


//    int num=atoi(argv[1]);
    if(argc<5){
        cout<<"Usage: "<<argv[0]<<"<rmap_file_in_val_format> <cluster_centers_file> <total_streams> <stream_index>"<<endl;
        return 1;
    }

    std::ifstream outputfile(argv[2]); // cluster centers

    int tot_streams=atoi(argv[3]); // total number of streams
    int index=atoi(argv[4]); // index of current stream

    string str;
    time_t ostart = clock();
    string sttt;
    int maxcluster=0;
    int cl;

    std::vector< std::vector<float> > readmeans;
    int cntmeans=0;

    while(cntmeans<4){
        getline(outputfile,str);
        std::stringstream sss(str);
        std::vector<float> temp;
        float temfl;

        while(sss>>temfl)
            temp.push_back(temfl);

        readmeans.push_back(temp);
        cntmeans++;
    }



    read_data(argv[1]); // rmap file

    vector<int> cluster_assignment;

    pointVec points;

    for(int i=0;i<readmeans[0].size();i++){

        point_t newpt;
        for(int ii=0;ii<readmeans.size();ii++){
            newpt.push_back(readmeans[ii][i]*1000);
        }
        points.push_back(newpt);
    }


    maxcluster = points.size();
    KDTree tree(points);

    for(int i=0;i<data.size();i++){

        point_t newpt;
        for(int ii=0;ii<data[i].size();ii++){
            newpt.push_back(data[i][ii]*1000);
        }
        cluster_assignment.push_back(tree.nearest_index_nod2(newpt,newpt));

        if(i<10)
            cout<<cluster_assignment.back()<<" "<<rmap_belong[i]<<endl;
     }


    int rmap_cnt=rmap_count;

    int eachstream = (rmap_count/tot_streams);

    int start=eachstream*(index-1);
    int fin=eachstream*(index);

    cout<<rmap_cnt<<" "<<tot_streams<<" "<<index<<endl;

    cout<<"\n Rmap start index:"<< start<<" Rmap end index:"<<fin<<"\n";


    vector < vector <int> > rmap_to_cluster;

    vector < vector <int> > cluster_to_rmap;

    vector<int> temp;

    cluster_to_rmap.resize(maxcluster,temp);

    int lineno=0;

    rmap_to_cluster.resize(rmap.size(),temp);

    for(int i=0;i<cluster_assignment.size();i++){

        if(cluster_assignment[i]<0 || cluster_assignment[i]>=cluster_to_rmap.size())
            cout<<"error here "<<" "<<cluster_to_rmap.size()<<endl;

        cluster_to_rmap[cluster_assignment[i]].push_back(nameid.find(rmap_belong[i])->second);
        rmap_to_cluster[nameid.find(rmap_belong[i])->second].push_back(cluster_assignment[i]);
    }


	time_t reading_t_end = clock();
	int numrmaps = rmap.size();

	cout<<"Data read time: "<<double(reading_t_end - ostart) / CLOCKS_PER_SEC << endl;

    cout<<"Now finding alignments with related Rmaps "<<" "<< rmap.size()<< " "<<rmap[0].size()<<endl;

    ofstream numrels("numrels");

    int dval=4;

    for(int i=start;i<min(rmap_cnt,fin);i++){

        vector<int> relation_vector(numrmaps,0);
        vector<int> collect_relations;
        int cntrels=0;

        for(int j=0;j<rmap_to_cluster[i].size();j++){

            int cluster=rmap_to_cluster[i][j];
            for(int k=0;k<cluster_to_rmap[cluster].size();k++){

                int otherrmap = cluster_to_rmap[cluster][k];
                relation_vector[otherrmap]++;

                if(relation_vector[otherrmap]==dval && otherrmap!=i){
                    cntrels++;
                     alignment = optimized_overlap_alignment(rmap[i],rmap[otherrmap],0,0);

                    if(alignment[0].first>8 && alignment[0].second>0){
                        rmap_relations_final[i].push_back(otherrmap);
                        rmap_alignments[i].push_back(alignment);
                    }

                    if(rmap_alignments[i].size()==max_alignments)
                        break;

                }


            }

            if(rmap_alignments[i].size()==max_alignments)
                        break;
        }

        multi_align_grid[i].resize(rmap_alignments[i].size());

        for(int k=0; k<multi_align_grid[i].size(); k++){
                multi_align_grid[i][k].resize(rmap[i].size(),make_trip(0,0,0));

        }

        cout<<i<<" "<<cntrels<<" "<<rmap_alignments[i].size()<<" "<<multi_align_grid[i].size()<<endl;

//        if(i%100==0)
            //cout<<i<<" "<<rmap_alignments[i].size()<<endl;

    }

    cout <<endl<< "Multi alignment grid building begins....."<<endl;

    build_multi_align_grid();

    time_t mag_end = clock();

//    print_multi_align_grid();
    cout << "Finding consensus....."<<endl;
    find_concensus(argv[4]);

    time_t find_con_end = clock();

    cout<<endl<<"Rmap correction begins.....";
//    rmaps_to_corr();
    correct_rmaps(argv[4]);

    time_t corr_rmaps_end = clock();

    time_t oend = clock();
    cout << endl << "total running time:	" << double(oend - ostart) / CLOCKS_PER_SEC << endl;

}


