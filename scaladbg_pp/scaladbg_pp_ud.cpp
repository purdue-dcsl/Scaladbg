/**
 * @file scaladbg_pp_ud.cpp : replace the idba_ud.cpp in idba/src/release/idba_ud.cpp
 * After doing the 'make' command, idba_ud is the resulting executable file
 * ScalaDBG using idba_ud.cpp as a base for starting implementation. 
 * based on Yu Peng idba_ud.cpp implementation
 */

#include <queue>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h> 
#define MASTER 0 
#define SEND 0
#define RECV 1
#define NOTHING 2 
#define INVALID -1
#define DEBUG 0
#define TIMING 1
#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>
#include "assembly/assembly_utility.h"
#include "assembly/local_assembler.h"
#include "basic/bit_operation.h"
#include "basic/histgram.h"
#include "graph/contig_graph.h"
#include "graph/hash_graph.h"
#include "graph/scaffold_graph.h"
#include "misc/hash_aligner.h"
#include "misc/log.h"
#include "misc/options_description.h"
#include "misc/utils.h"
#include "sequence/read_library.h"
#include "sequence/sequence.h"
#include "sequence/sequence_io.h"
#include "sequence/short_sequence.h"

using namespace std;

struct IDBAOption
{
    string directory;
    string read_file;
    string long_read_file;
    deque<string> extra_read_files;
    int mink;
    int maxk;
    int step;
    int inner_mink;
    int inner_step;
    int prefix_length;
    int min_count;
    int min_support;
    int min_contig;
    double similar;
    int max_mismatch;
    int seed_kmer_size;
    int num_threads;
    int min_pairs;
    int max_gap;
    bool is_no_bubble;
    bool is_no_local;
    bool is_no_coverage;
    bool is_no_correct;
    bool is_pre_correction;
    string reference;

    IDBAOption()
    {
        extra_read_files.resize(4);
        directory = "out";
        mink = 20;
        maxk = 100;
        step = 20;
        inner_mink = 10;
        inner_step = 5;
        prefix_length = 3;
        min_count = 2;
        min_support = 1;
        min_contig = 200;
        similar = 0.95;
        max_mismatch = 3;
        seed_kmer_size = 30;
        num_threads = 0;
        min_pairs = 3;
        max_gap = 50;
        is_no_bubble = false;
        is_no_local = false;
        is_no_coverage = false;
        is_no_correct = false;
        is_pre_correction = false;
    }

    string log_file()
    { return directory + "/log"; }

    string kmer_file()
    { return directory + "/kmer"; }

    string align_file(int kmer_size)
    { return directory + FormatString("/align-%d", kmer_size); }

    string graph_file(int kmer_size)
    { return directory + FormatString("/graph-%d.fa", kmer_size); }

    string contig_file(int kmer_size)
    { return directory + FormatString("/contig-%d.fa", kmer_size); }

    string contig_info_file(int kmer_size)
    { return directory + FormatString("/contig-info-%d.fa", kmer_size); }

    string local_contig_file(int kmer_size)
    { return directory + FormatString("/local-contig-%d.fa", kmer_size); }

    string contig_file()
    { return directory + "/contig.fa"; }

    string scaffold_file(int level = 0)
    { return directory + (level == 0 ? "/scaffold.fa" : FormatString("/scaffold-level-%d.fa", level+1)); }

    string ref_contig_file()
    { return directory + "/ref_contig.fa"; }
};
AssemblyInfo assembly_info;
IDBAOption option;
double median = 0;
double sd = 0;
int read_length = 0;
static void newBuildHashGraph(int kmer_size, HashGraph& hash_graph);
static void UpdateHashGraph(int old_kmer_size,int kmer_size, HashGraph& hash_graph, int min_support, int rank);
static void get_max(std::vector<int> &ranks_left, std::vector<int> &curr_kvals, int &idx);
static void build_matrices(int size, int numkmers, std::vector<std::vector<int> >&send_recv,
    std::vector<std::vector<int> >&communication,std::vector<std::vector<int> >&k_val_work, 
    std::queue<int> kvals);
static void print_vectors(vector< vector <int> > &k_val_work,vector< vector <int> >&send_recv,
    vector< vector <int> >&communication);
void Assemble(HashGraph &hash_graph);
void AlignReads(const string &contig_file, const string &align_file);
void CorrectReads(int kmer_size);
void LocalAssembly(int kmer_size, int new_kmer_size, int rank);
void Scaffold(int kmer_size, int min_contig);
void AddPairs(int level, ScaffoldGraph &scaffold_graph, const string &read_file, const string &align_file);
void AlignReads(const string &contig_file, ShortReadLibrary &library, const string &align_file);

<<<<<<< HEAD:mpiBig/idba_ud.cpp
int main(int argc, char *argv[])
{
    double readTime, precorTime, recvTime;
    double   sendTime,totalTime;
    // double assemTime,alignTime, corrTime, scaffTime,buildTime,updateTime;
=======
int main(int argc, char *argv[]){
    double readTime, recvTime;
    double sendTime,totalTime;
>>>>>>> testing:scaladbg_pp/scaladbg_pp_ud.cpp
    int rank;
    int size;
    OptionsDescription desc;
    
    desc.AddOption("out", "o", option.directory, "output directory");
    desc.AddOption("read", "r", option.read_file, FormatString("fasta read file (<=%d)", ShortSequence::max_size()));
    desc.AddOption("read_level_2", "", option.extra_read_files[0], "paired-end reads fasta for second level scaffolds");
    desc.AddOption("read_level_3", "", option.extra_read_files[1], "paired-end reads fasta for third level scaffolds");
    desc.AddOption("read_level_4", "", option.extra_read_files[2], "paired-end reads fasta for fourth level scaffolds");
    desc.AddOption("read_level_5", "", option.extra_read_files[3], "paired-end reads fasta for fifth level scaffolds");
    desc.AddOption("long_read", "l", option.long_read_file, FormatString("fasta long read file (>%d)", ShortSequence::max_size()));
    desc.AddOption("mink", "", option.mink, FormatString("minimum k value (<=%d)", Kmer::max_size()));
    desc.AddOption("maxk", "", option.maxk, FormatString("maximum k value (<=%d)", Kmer::max_size()));
    desc.AddOption("step", "", option.step, "increment of k-mer of each iteration");
    desc.AddOption("inner_mink", "", option.inner_mink, "inner minimum k value");
    desc.AddOption("inner_step", "", option.inner_step, "inner increment of k-mer");
    desc.AddOption("prefix", "", option.prefix_length, "prefix length used to build sub k-mer table");
    desc.AddOption("min_count", "", option.min_count, "minimum multiplicity for filtering k-mer when building the graph");
    desc.AddOption("min_support", "", option.min_support, "minimum supoort in each iteration");
    desc.AddOption("num_threads", "", option.num_threads, "number of threads");
    desc.AddOption("seed_kmer", "", option.seed_kmer_size, "seed kmer size for alignment");
    desc.AddOption("min_contig", "", option.min_contig, "minimum size of contig");
    desc.AddOption("similar", "", option.similar, "similarity for alignment");
    desc.AddOption("max_mismatch", "", option.max_mismatch, "max mismatch of error correction");
    desc.AddOption("min_pairs", "", option.min_pairs, "minimum number of pairs");
    desc.AddOption("no_bubble", "", option.is_no_bubble, "do not merge bubble");
    desc.AddOption("no_local", "", option.is_no_local, "do not use local assembly");
    desc.AddOption("no_coverage", "", option.is_no_coverage, "do not iterate on coverage");
    desc.AddOption("no_correct", "", option.is_no_correct, "do not do correction");
    desc.AddOption("pre_correction", "", option.is_pre_correction, "perform pre-correction before assembly");
    
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if(DEBUG)printf("My rank is %i and size is %i\n",rank,size);

    try{
        desc.Parse(argc, argv);

        if (option.read_file == "" && option.long_read_file == "")
            throw logic_error("not enough parameters");

        if (option.maxk < option.mink)
            throw invalid_argument("mink is larger than maxk");

        if (option.maxk > (int)Kmer::max_size())
            throw invalid_argument("maxk is too large");
    }
    catch (exception &e){
        cerr << e.what() << endl;
        cerr << "ScalaDBG_PP" << endl;
        cerr << "Usage: idba_ud -r read.fa -o output_dir" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }
    if(rank==MASTER)
        MakeDir(option.directory);
    MPI_Barrier(MPI_COMM_WORLD);//so we can get just the computation time, not setup.
    totalTime = -MPI_Wtime();

    LogThread log_thread(option.log_file());

    string begin_file = option.directory + "/begin";
    fclose(OpenFile(begin_file, "wb"));

    if (option.num_threads == 0)
        option.num_threads = omp_get_max_threads();
    else
        omp_set_num_threads(option.num_threads);
    
    if(rank==MASTER)
        printf("number of threads %d\n",option.num_threads);
    
    readTime = -MPI_Wtime();
    ReadInput(option.read_file, option.long_read_file, assembly_info);
    deque<Sequence> extra_reads;
    for (unsigned i = 0; i < option.extra_read_files.size(); ++i){
        if (option.extra_read_files[i] != ""){
            deque<Sequence> reads;
            ReadSequence(option.extra_read_files[i], reads);
            extra_reads.insert(extra_reads.end(), reads.begin(), reads.end());
        }
    }
    readTime += MPI_Wtime();
    if(rank==MASTER){
        cout << "reads " << assembly_info.reads.size() << endl;
        cout << "long reads " << assembly_info.long_reads.size() << endl;
        cout << "extra reads " << extra_reads.size() << endl;
        if(TIMING)printf("Read Time (ReadSequence) was %f seconds\n",readTime);
    }
    assembly_info.long_reads.insert(assembly_info.long_reads.end(), extra_reads.begin(), extra_reads.end());
    assembly_info.ClearStatus();

    read_length = assembly_info.read_length();
    if(rank==MASTER)
        cout << "read_length " << read_length << endl;

    vector< vector <int> >k_val_work; 
    vector< vector <int> > send_recv; 
    vector<vector<int> >communication;
    int numkmers = (option.maxk - option.mink)/option.step + 1;
    build_matrices(size, numkmers, send_recv, communication, k_val_work);

    if(DEBUG==1 && rank == MASTER)
        print_vectors(k_val_work,send_recv,communication);
    
    //this is the start of the computation loop
    int kmer_size;
    for(unsigned int i = 0; i < send_recv[rank].size(); i++){
        if(send_recv[rank][i] == SEND){
            if(communication[rank][i] != -1){
                int sendingTo = communication[rank][i];//should reduce
                sendTime = -MPI_Wtime();
                MPI_Send(&kmer_size,1,MPI_INT,sendingTo,0,MPI_COMM_WORLD);
                sendTime += MPI_Wtime();
                if(TIMING)printf("Sending from %i to %i took %f seconds\n",rank, sendingTo, sendTime);               
            }
            if(!k_val_work[rank].empty()){
                //build next value, if in here we have more work
                double build_time = -MPI_Wtime();
                kmer_size = option.mink + k_val_work[rank][0]*option.step;
                k_val_work[rank].erase(k_val_work[rank].begin());//pop off k_value we are going to build
                HashGraph hash_graph(kmer_size);
                newBuildHashGraph(kmer_size,hash_graph);
<<<<<<< HEAD:mpiBig/idba_ud.cpp
                newbuildTime +=MPI_Wtime();
                printf("New Build time for rank %i nad ksize %i took %f seconds\n",rank, kmer_size, newbuildTime);
=======
                
>>>>>>> testing:scaladbg_pp/scaladbg_pp_ud.cpp
                if((i+1) < send_recv[rank].size() && send_recv[rank][i+1] == RECV){//if we are receiving next we will patch, need contigs from this one to patch :)
                    Assemble(hash_graph);
                }
                else{
                    char filename [33];
                    ofstream outfile;
                    sprintf (filename,"%s/%d",option.directory.c_str(),kmer_size);
                    outfile.open(filename, ios::binary | ios::out);
                    outfile << hash_graph;
                    outfile.close();
<<<<<<< HEAD:mpiBig/idba_ud.cpp
                    writingTime += MPI_Wtime();
                    printf("Write time for rank %i and ksize %i took %f seconds\n",rank, kmer_size, writingTime);
=======
>>>>>>> testing:scaladbg_pp/scaladbg_pp_ud.cpp
                }
                build_time += MPI_Wtime();
                if(TIMING)printf("Build time for kmer_size %i is %f seconds\n",kmer_size,build_time);
            }
        }
        else if(send_recv[rank][i] == RECV){
            int rec_kmer_size;
            recvTime = -MPI_Wtime();
            int source = communication[rank][i];
            MPI_Recv(&rec_kmer_size,1,MPI_INT,source,0, MPI_COMM_WORLD, &status);
            recvTime += MPI_Wtime();
<<<<<<< HEAD:mpiBig/idba_ud.cpp
            printf("Rank %i Receiving from %i took %f seconds, kmer_size recv is %i\n",rank, source,recvTime, rec_kmer_size);

            // double alignTime = -MPI_Wtime();
            // AlignReads(option.contig_file(kmer_size), option.align_file(kmer_size));
            // alignTime += MPI_Wtime();
            // double corrTime = -MPI_Wtime();
            // CorrectReads(kmer_size);
            // corrTime += MPI_Wtime();
            // printf("Align Reads for kmer size %i took %f seconds\n", kmer_size,alignTime);
            // printf("Correct Reads for kmer size %i took %f seconds\n", kmer_size,corrTime);
            // assembly_info.ClearStatus();
            // double local_time = -MPI_Wtime();
            // LocalAssembly(kmer_size, rec_kmer_size);
            // local_time +=MPI_Wtime();
            // printf("Rank %i Local Assembly with new ksize %i took %f seconds\n",rank, rec_kmer_size,local_time);
            double patch_time = -MPI_Wtime();

=======
            if(TIMING)printf("Rank %i Receiving from %i took %f seconds, kmer_size recv is %i\n",rank, source,recvTime, rec_kmer_size);
            double align_time = -MPI_Wtime();
            AlignReads(option.contig_file(kmer_size), option.align_file(kmer_size));
            align_time += MPI_Wtime();
            if(TIMING)printf("AlignReads for kmer_size %i took %f seconds\n",kmer_size,align_time);
            double correct_time = -MPI_Wtime();
            CorrectReads(kmer_size);
            correct_time+=MPI_Wtime();
            if(TIMING)printf("CorrectReads for kmer_size %i took %f seconds\n",kmer_size,correct_time);
            assembly_info.ClearStatus();
            
            double local_time = -MPI_Wtime();
            LocalAssembly(kmer_size, rec_kmer_size,rank);
            local_time+=MPI_Wtime();
            if(TIMING)printf("LocalAssembly for kmer_size %i took %f seconds\n",kmer_size,local_time);


            double patch_time = -MPI_Wtime();
>>>>>>> testing:scaladbg_pp/scaladbg_pp_ud.cpp
            //patch now
            char filename [33];
            int old_kmer_size = kmer_size;
            kmer_size = rec_kmer_size;
            HashGraph new_hash_graph(kmer_size);
            ifstream infile;
<<<<<<< HEAD:mpiBig/idba_ud.cpp
            sprintf (filename,"%d",kmer_size);
            infile.open(filename, ios::binary | ios::in);
            infile>>new_hash_graph;
            infile.close();
            
            UpdateHashGraph(old_kmer_size,kmer_size,new_hash_graph,option.min_count);
=======
            sprintf (filename,"%s/%d",option.directory.c_str(),kmer_size);
            infile.open(filename, ios::binary | ios::in);
            infile>>new_hash_graph;
            infile.close();
            UpdateHashGraph(old_kmer_size,kmer_size,new_hash_graph,option.min_count, rank);
>>>>>>> testing:scaladbg_pp/scaladbg_pp_ud.cpp
            patch_time += MPI_Wtime();
            if(TIMING)printf("Patching or UpdateHashGraph for kmer_size %i took %f seconds\n",kmer_size,patch_time);
            if(((i+1) < send_recv[rank].size() && send_recv[rank][i+1] == RECV)){//if we are receiving next we will patch, need contigs from this one to patch :)
                double assemble_time = -MPI_Wtime();
                Assemble(new_hash_graph);
                assemble_time += MPI_Wtime();
                if(TIMING)printf("Assemble for kmer_size %i took %f seconds\n",kmer_size,assemble_time);   
            }

            if(rank==MASTER && (i+1)==send_recv[rank].size()){
                double assemble_time = -MPI_Wtime();
                Assemble(new_hash_graph);
                assemble_time += MPI_Wtime();
                if(TIMING)printf("Assemble for kmer_size %i took %f seconds\n",kmer_size,assemble_time);
            }
        }
        else{
            ;
        }
    }
    if(rank==MASTER){
        if(DEBUG)printf("rank %i is the master finishing, kmer_size is %i, option.maxk is %i\n",rank,kmer_size,option.maxk);
        readTime = -MPI_Wtime();
        kmer_size = option.maxk;
        deque<Sequence> contigs;
        deque<string> names;
        ReadSequence(option.contig_file(kmer_size), contigs, names);
        FastaWriter writer(option.contig_file());
        for (unsigned i = 0; i < contigs.size(); ++i){
            if ((int)contigs[i].size() >= option.min_contig)
                writer.Write(contigs[i], names[i]);
        }
        readTime += MPI_Wtime();
        if(TIMING)printf("ReadSequence/writer time for kmer_size %i took %f seconds\n",kmer_size,readTime);
        // scaffTime = -MPI_Wtime();
        // Scaffold(option.maxk, option.min_contig);
        // scaffTime += MPI_Wtime();
        // printf("Scaffold took %f seconds\n", scaffTime);
        string end_file = option.directory + "/end";
        fclose(OpenFile(end_file, "wb"));
    }
    totalTime += MPI_Wtime();
    if(rank==MASTER)
        printf("Total Program time was %f seconds\n", totalTime);
    fflush(stdout);
    MPI_Finalize();
    return 0;
}

static void print_vectors(vector< vector <int> > &k_val_work,vector< vector <int> >&send_recv,
    vector< vector <int> >&communication){
    printf("k_val_work\n");
    for(unsigned int i = 0; i < k_val_work.size(); i++){
        for(unsigned int j =0; j < k_val_work[i].size();j++){
            printf("%i ", k_val_work[i][j]);;
        }
        printf("\n");
    }

    printf("send_recv\n");
    for(unsigned int i = 0; i < send_recv.size(); i++){
        for(unsigned int j =0; j < send_recv[i].size();j++){
            printf("%i ", send_recv[i][j]);;
        }
        printf("\n");
    }

    printf("communication\n");
    for(unsigned int i = 0; i < communication.size(); i++){
        for(unsigned int j =0; j < communication[i].size();j++){
            printf("%i ", communication[i][j]);;
        }
        printf("\n");
    }
}

static void get_max(std::vector<int> &ranks_left, std::vector<int> &curr_kvals, int &idx){
    int max_ele = 0;
    for (unsigned int i = 0; i < ranks_left.size(); i++){
        if(ranks_left[i]>=max_ele){
            if(ranks_left[i]>max_ele){
                max_ele = ranks_left[i];
                idx = i;
            }
            else if(ranks_left[i]==max_ele && curr_kvals[i]<curr_kvals[idx]){
                max_ele = ranks_left[i];
                idx = i;
            }
            else{
                ;
            }
        }
    }
    if(idx!=INVALID)
        ranks_left[idx] = INVALID;
}

static void build_matrices(int size, int numkmers, std::vector<std::vector<int> >&send_recv,std::vector<std::vector<int> >&communication,std::vector<std::vector<int> >&k_val_work){
    vector<int> rank_weight(size,0);//initialize to number of mpi_processes, all with 0 weight
    vector<int> curr_kvals(size,0);
    vector<int> topush;
    queue<int> kvals;
    for(int i = 0; i < numkmers; i++){
        kvals.push(i);
    }
    for(int i = 0; i < size; i++){
        curr_kvals[i]=-1;
        k_val_work.push_back(topush);
        send_recv.push_back(topush);
        communication.push_back(topush);
    }
    //get the first round initialization, at least 1 round always
    for(int i = 0; i < size; i++){
        if(!kvals.empty()){
            k_val_work[i].push_back(kvals.front());
            kvals.pop();
            rank_weight[i]=1;
            send_recv[i].push_back(0);
        }
        else{
            send_recv[i].push_back(2);
        }
        communication[i].push_back(-1);
    }
    while(!kvals.empty() || rank_weight[0]!=numkmers){
        //do a round
        vector<int> ranks_left(rank_weight.begin(), rank_weight.end());
        vector<bool> visited(size,false);
        for(int i = 0; i < size/2; i++){
            //get the process with the most kvalues processed.
            int f_idx=-1,s_idx=-1;
            get_max(ranks_left,curr_kvals,f_idx);
            //get the process with the second most kvalues processed.
            get_max(ranks_left,curr_kvals,s_idx);
            //the one with the most kvalues receives from the one with the second most.
            if(f_idx != INVALID && s_idx != INVALID){
                visited[f_idx]=true;
                visited[s_idx]=true;
                if(rank_weight[f_idx]>0 && rank_weight[s_idx]>0){//both have been working on something
                    rank_weight[f_idx]+=rank_weight[s_idx];
                    rank_weight[s_idx]=-1;
                    communication[f_idx].push_back(s_idx);//this will receive
                    communication[s_idx].push_back(f_idx);//this will send
                    send_recv[f_idx].push_back(1);//this will receive
                    send_recv[s_idx].push_back(0);//this will send
                    curr_kvals[f_idx]=curr_kvals[s_idx];
                    curr_kvals[s_idx] = -1;
                }
            }
            else{
                for(int j = 0; j < size; j++){
                    if(visited[j]==false){
                       communication[j].push_back(-1);
                       send_recv[j].push_back(2); 
                       visited[j]=true;
                    }
                }
                //need to fill in the last entries in send_recv and communication with do nothing and -1
            }
        }
        for(int i = 0; i < size; i ++){
            if(send_recv[i].back()==0 && !kvals.empty()){
                curr_kvals[i]=kvals.front();
                k_val_work[i].push_back(kvals.front());
                kvals.pop();
                rank_weight[i]=1;
            }
        }
    }
}


static void newBuildHashGraph(int kmer_size, HashGraph& hash_graph){
    char buffer [33];
    sprintf (buffer,"%d",kmer_size);
    BuildKmerFile(assembly_info, kmer_size, option.min_count, option.prefix_length, option.kmer_file()+buffer);
    ReadKmerFile(option.kmer_file()+buffer, hash_graph);    
    hash_graph.RefreshEdges();
    InsertInternalKmers(assembly_info, hash_graph, option.min_count);
}

static void UpdateHashGraph(int old_kmer_size,int kmer_size, HashGraph& hash_graph, int min_count, int rank){
    deque<Sequence> contigs;
    ReadSequence(option.contig_file(old_kmer_size), contigs);
    
    deque<Sequence> local_contigs;//used only when LocalAssembly is used
    ReadSequence(option.local_contig_file(old_kmer_size), local_contigs);//used only when LocalAssembly is used
    
    deque<Sequence> multi_contigs;
    ReadSequence(option.graph_file(old_kmer_size), multi_contigs);
    deque<Sequence> old_contigs;
    old_contigs.insert(old_contigs.end(), contigs.begin(), contigs.end());
    old_contigs.insert(old_contigs.end(), local_contigs.begin(), local_contigs.end());//used only when LocalAssembly is used
    old_contigs.insert(old_contigs.end(), multi_contigs.begin(), multi_contigs.end());
    contigs.clear();
    local_contigs.clear();//used only when LocalAssembly is used
    multi_contigs.clear();
    
    #pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)old_contigs.size(); ++i)
        hash_graph.InsertUncountKmers(old_contigs[i]);
    hash_graph.RestoreAndMergeEdges();
    hash_graph.RefreshEdges();
}

void Assemble(HashGraph &hash_graph){
    cout << "kmers " << hash_graph.num_vertices() << " "<< hash_graph.num_edges() << endl;

    int kmer_size = hash_graph.kmer_size();
    double min_cover = max(1, (kmer_size == option.mink ? option.min_count : option.min_support));

    Histgram<int> hist = hash_graph.coverage_histgram();
    double expected_coverage = hist.mean();

    deque<Sequence> contigs;
    deque<ContigInfo> contig_infos;
    hash_graph.Assemble(contigs, contig_infos);
    hash_graph.clear();

    {
        HashGraph tmp_hash_graph;
        tmp_hash_graph.swap(hash_graph);
    }

    ContigGraph contig_graph(kmer_size, contigs, contig_infos);
    contigs.clear();
    contig_infos.clear();

    contig_graph.RemoveDeadEnd(option.min_contig);

    if (!option.is_no_bubble) 
    {
        int bubble = contig_graph.RemoveBubble();
        cout << "merge bubble " << bubble << endl;
        contig_graph.MergeSimilarPath();
    }

    if (!option.is_no_coverage)
        contig_graph.RemoveLocalLowCoverage(min_cover, option.min_contig, 0.1);

    contig_graph.SortVertices();
    contig_graph.GetContigs(contigs, contig_infos);
    WriteSequence(option.graph_file(kmer_size), contigs);
    contigs.clear();
    contig_infos.clear();

    if (!option.is_no_coverage)
    {
        double ratio = (kmer_size < option.maxk) ? 0.5 : 0.2;
        if (ratio < 2.0 / expected_coverage)
            ratio = 2.0 / expected_coverage;
        contig_graph.IterateLocalCoverage(option.min_contig, ratio, min_cover, 1e100, 1.1);
        contig_graph.MergeSimilarPath();
    }

    deque<Sequence> multi_contigs;
    deque<ContigInfo> multi_contig_infos;
    contig_graph.SortVertices();
    contig_graph.GetContigs(multi_contigs, multi_contig_infos);
    PrintN50(multi_contigs);
    WriteContig(option.contig_file(kmer_size), multi_contigs, multi_contig_infos, FormatString("contig-%d", kmer_size));
}

void AlignReads(const string &contig_file, const string &align_file){
    deque<Sequence> contigs;
    ReadSequence(contig_file, contigs);

    HashAligner hash_aligner(option.seed_kmer_size, option.min_contig, 2);
    hash_aligner.Initialize(contigs);

    int64_t num_aligned_reads = AlignReads(assembly_info, hash_aligner, option.similar, align_file, true);
    cout << "aligned " << num_aligned_reads << " reads" << endl;
}

void CorrectReads(int kmer_size){
    if (option.is_no_correct)
        return;

    deque<Sequence> contigs;
    deque<string> names;
    deque<ContigInfo> contig_infos;
    ReadSequence(option.contig_file(kmer_size), contigs, names);
    CorrectReads(assembly_info, contigs, contig_infos, option.align_file(kmer_size), option.max_mismatch);
    WriteContig(option.contig_file(kmer_size), contigs, contig_infos, FormatString("contig-%d", kmer_size));
}

void LocalAssembly(int kmer_size, int new_kmer_size, int rank){
    printf("Rank %i at line %i\n",rank,__LINE__);
    printf("In LocalAssembly with kmer_size %i and new_kmer_size %i\n",kmer_size,new_kmer_size);
    EstimateDistance(option.align_file(kmer_size), median, sd);
    if (median < 0 || median != median || sd != sd || sd > 2*median){
        cout << "invalid insert distance" << endl;
        deque<Sequence> local_contigs;
        WriteSequence(option.local_contig_file(kmer_size), local_contigs, FormatString("local_contig_%d", kmer_size));
        return;
    }

    deque<ShortSequence> &reads = assembly_info.reads;

    deque<Sequence> contigs;
    ReadSequence(option.contig_file(kmer_size), contigs);

    LocalAssembler local_assembler;
    local_assembler.Initialize(assembly_info, contigs);
    local_assembler.set_num_threads(option.num_threads);
    local_assembler.set_mink(option.inner_mink);
    local_assembler.set_maxk(new_kmer_size);
    local_assembler.set_step(option.inner_step);
    local_assembler.set_min_contig(option.min_contig);
    local_assembler.set_insert_distance(median, sd);

    FILE *falign = OpenFile(option.align_file(kmer_size), "rb");
    int buffer_size = (1 << 20) * option.num_threads;
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size){
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);

        ReadHashAlignerRecordBlock(falign, all_records);
#pragma omp parallel for
        for (int i = 0; i < size; ++i){
            HashAlignerRecord &record = all_records[i];

            if (record.match_length != 0)
                local_assembler.AddReadByHashAlignerRecord(record, offset + i);
        }
    }
    fclose(falign);

    deque<Sequence> local_contigs;

    if (!option.is_no_local)
        local_assembler.Assemble(local_contigs);
    
    int num_seed_contigs = 0;
    for (unsigned i = 0; i < contigs.size(); ++i){
        if ((int)contigs[i].size() > option.min_contig)
            ++num_seed_contigs;
    }

    cout << "seed contigs " << num_seed_contigs <<  " local contigs " << local_contigs.size() << endl;
    WriteSequence(option.local_contig_file(kmer_size), local_contigs, FormatString("local_contig_%d", kmer_size));
}

void Scaffold(int kmer_size, int min_contig){
    assembly_info.reads.clear();
    assembly_info.long_reads.clear();

    deque<Sequence> contigs;
    ReadSequence(option.contig_file(option.maxk), contigs);
    ScaffoldGraph scaffold_graph(option.maxk, contigs);

    deque<string> read_files;
    read_files.push_back(option.read_file);
    for (unsigned i = 0; i < option.extra_read_files.size(); ++i){
        if (option.extra_read_files[i] != "")
            read_files.push_back(option.extra_read_files[i]);
    }

    for (int level = 0; level < (int)read_files.size(); ++level)
        AddPairs(level, scaffold_graph, read_files[level], option.align_file(option.maxk) + FormatString("-%d", level));

    for (int level = 0; level < (int)read_files.size(); ++level){
        scaffold_graph.BuildEdges();
        scaffold_graph.FilterEdges(option.min_pairs, scaffold_graph.sd(level) * 4);
        scaffold_graph.ParseEdges();

        cout << "edgs " << scaffold_graph.num_edges(level) << endl;
        scaffold_graph.RemoveTransitiveConnections(level);

        deque<ContigGraphPath> paths;
        scaffold_graph.Assemble(level, paths);

        deque<Sequence> contigs;
        scaffold_graph.Assemble(level, contigs);
        PrintN50(contigs);

        WriteSequence(option.scaffold_file(level), contigs, "scaffold");

        scaffold_graph.Initialize(paths);
    }
}

void AddPairs(int level, ScaffoldGraph &scaffold_graph, const string &read_file, const string &align_file){
    ShortReadLibrary short_read_library;
    ReadLibrary(read_file, short_read_library);
    cout << "reads " << short_read_library.size() << endl;
    AlignReads(option.contig_file(option.maxk), short_read_library, align_file);

    EstimateDistance(align_file, median, sd);

    deque<Sequence> contigs;
    ReadSequence(option.contig_file(option.maxk), contigs);

    deque<ContigInfo> contig_infos(contigs.size());
    vector<int> num_aligned_reads(contigs.size(), 0);
    vector<double> coverage(contigs.size());

    deque<ShortSequence> &reads = short_read_library.reads();

    FILE *falign = OpenFile(align_file, "rb");
    int buffer_size = (1 << 20) * option.num_threads;
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size){
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);

        ReadHashAlignerRecordBlock(falign, all_records);
#pragma omp parallel for
        for (int i = 0; i < size; ++i){
            HashAlignerRecord &record = all_records[i];

            if (record.match_length != 0){
#pragma omp atomic
                ++num_aligned_reads[record.ref_id];
            }
        }
    }
    fclose(falign);

    double sum_coverage = 0;
    double sum_length = 0;
#pragma omp parallel for reduction(+: sum_coverage, sum_length)
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i){
        if ((int)contigs[i].size() > option.min_contig)
        {
            sum_coverage += num_aligned_reads[i];
            sum_length += contigs[i].size() - read_length + 1;
            coverage[i] = 1.0 * num_aligned_reads[i] / (contigs[i].size() - reads[0].size() + 1);
            contig_infos[i].set_kmer_count(num_aligned_reads[i]);
        }
    }
    double mean_coverage = sum_coverage / sum_length;
    cout << "expected coverage " << mean_coverage << endl;

    int num_connections = 0;
    falign = OpenFile(align_file, "rb");
    for (unsigned i = 0; i < reads.size(); i += 2){
        deque<HashAlignerRecord> records1;
        deque<HashAlignerRecord> records2;
        ReadHashAlignerRecords(falign, records1);
        ReadHashAlignerRecords(falign, records2);

        for (unsigned j = 0; j < records1.size(); ++j){
            for (unsigned k = 0; k < records2.size(); ++k){
                HashAlignerRecord &r1 = records1[j];
                HashAlignerRecord &r2 = records2[k];
                r2.ReverseComplement();

                if (r1.ref_length > option.min_contig && r2.ref_length > option.min_contig
                        && r1.ref_from - r1.query_from > r1.ref_length - median - 3*sd
                        && r2.ref_to + r2.query_length - r2.query_to < median + 3*sd
                        && r1.ref_id != r2.ref_id
                        ){
                    int d = median - (r1.ref_length - (r1.ref_from - r1.query_from)) - (r2.ref_to + r2.query_length - r2.query_to);
                    scaffold_graph.AddPair(level, (r1.ref_id*2 + r1.is_reverse), (r2.ref_id*2 + r2.is_reverse), d);
                    ++num_connections;
                }
            }
        }
    }

    scaffold_graph.set_library_info(level, read_length, mean_coverage, median, sd);
}

void AlignReads(const string &contig_file, ShortReadLibrary &library, const string &align_file){
    deque<Sequence> contigs;
    ReadSequence(contig_file, contigs);

    assembly_info.reads.swap(library.reads());
    HashAligner hash_aligner(option.seed_kmer_size, option.min_contig, 2);
    hash_aligner.Initialize(contigs);

    int64_t num_aligned_reads = AlignReads(assembly_info, hash_aligner, option.similar, align_file, true);
    cout << "aligned " << num_aligned_reads << " reads" << endl;

    assembly_info.reads.swap(library.reads());
}
