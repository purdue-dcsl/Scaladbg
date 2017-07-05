/**
 * @file scaladbg_sp_ud.cpp : replace the idba_ud.cpp in idba/src/release/idba_ud.cpp
 * After doing the 'make' command, idba_ud is the resulting executable file
 * ScalaDBG using idba_ud.cpp as a base for starting implementation. 
 * based on Yu Peng idba_ud.cpp implementation
 */

#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <mpi.h> 
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

#define MASTER 0
#define DEBUG 0
#define TIMING 0

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
static void UpdateHashGraph(int old_kmer_size,int kmer_size, HashGraph& hash_graph, int min_support);
void Assemble(HashGraph &hash_graph);
void AlignReads(const string &contig_file, const string &align_file);
void CorrectReads(int kmer_size);
void LocalAssembly(int kmer_size, int new_kmer_size);
void Scaffold(int kmer_size, int min_contig);
void AddPairs(int level, ScaffoldGraph &scaffold_graph, const string &read_file, const string &align_file);
void AlignReads(const string &contig_file, ShortReadLibrary &library, const string &align_file);

int main(int argc, char *argv[])
{
    double totalTime;
    int rank,size;
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
    
    if(size<2){
        printf("Only have one node, run with original IDBA_UD\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    
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
        cerr << "ScalaDBG_SP" << endl;
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
    if(rank==MASTER && DEBUG)
        printf("number of threads %d\n",option.num_threads);

    double readTime = -MPI_Wtime();

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
        printf("Read Time (ReadSequence) was %f seconds\n",readTime);
    }

    assembly_info.long_reads.insert(assembly_info.long_reads.end(), extra_reads.begin(), extra_reads.end());
    assembly_info.ClearStatus();

    read_length = assembly_info.read_length();
    if(rank==MASTER)
        cout << "read_length " << read_length << endl;

//This is the start of the computation loop
    //get the kmer_size I should compute
    int numkmers = (option.maxk - option.mink)/option.step + 1;
    if(numkmers != size){
        printf("This version doesn't work unless number of kmers equals number of nodes\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    int kmer_size = rank*option.step+option.mink;
    if(kmer_size<=option.maxk){
        double buildTime = -MPI_Wtime();
        HashGraph hash_graph(kmer_size);
        newBuildHashGraph(kmer_size,hash_graph);
        buildTime += MPI_Wtime();
        if(TIMING)printf("newBuildHashGraph for rank %i took %f seconds\n",rank,buildTime);
        if (rank==MASTER){
            int recv_rank;
            double assemTime = -MPI_Wtime();                   
            Assemble(hash_graph);
            assemTime += MPI_Wtime();
            if(TIMING)printf("MASTER Took %f seconds to Assemble kmer %i\n",assemTime,kmer_size);
            for(int i = 1; i < numkmers; i++){
                double recvTime = -MPI_Wtime();
                MPI_Recv(&recv_rank,1,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
                recvTime+=MPI_Wtime();
                if(TIMING)printf("MASTER recv from rank %i took %f seconds\n",recv_rank,recvTime);
            }
            //now we have received everything, continue!
            while (true){
                int old_kmer_size = kmer_size;
                kmer_size = min(option.maxk, kmer_size + option.step);
                if (old_kmer_size < option.maxk)
                {
                    double alignTime = -MPI_Wtime();
                    AlignReads(option.contig_file(old_kmer_size), option.align_file(old_kmer_size));
                    alignTime += MPI_Wtime();
                    double corrTime = -MPI_Wtime();
                    CorrectReads(old_kmer_size);
                    corrTime += MPI_Wtime();
                    if(TIMING)printf("Align Reads for kmer size %i took %f seconds\n", old_kmer_size,alignTime);
                    if(TIMING)printf("Correct Reads for kmer size %i took %f seconds\n", old_kmer_size,corrTime);
                    assembly_info.ClearStatus();
                    double local_time = -MPI_Wtime();
                    LocalAssembly(old_kmer_size, kmer_size);
                    local_time += MPI_Wtime();
                    if(TIMING)printf("LocalAssembly for old_kmer_size %i and kmer_size %i took %f seconds\n",
                        old_kmer_size,kmer_size,local_time);
                }
                double patch_time = -MPI_Wtime();
                HashGraph new_hash_graph(kmer_size);
                char filename [33];
                ifstream infile;
                sprintf (filename,"%s/%d",option.directory.c_str(), kmer_size);
                infile.open(filename, ios::binary | ios::in);
                infile>>new_hash_graph;
                infile.close();
                UpdateHashGraph(old_kmer_size,kmer_size,new_hash_graph,option.min_count);
                patch_time += MPI_Wtime();
                if(TIMING)printf("Patching or UpdateHashGraph for kmer_size %i took %f seconds\n",kmer_size,patch_time);
                double assemTime = -MPI_Wtime();                   
                Assemble(new_hash_graph);
                assemTime += MPI_Wtime();
                if(TIMING)printf("MASTER Took %f seconds to Assemble kmer %i\n",assemTime,kmer_size);
                if (kmer_size == option.maxk)
                    break;
            }
            double readTimehash_graph = -MPI_Wtime();
            int kmer_size = option.maxk;
            deque<Sequence> contigs;
            deque<string> names;
            ReadSequence(option.contig_file(kmer_size), contigs, names);
            FastaWriter writer(option.contig_file());
            for (unsigned i = 0; i < contigs.size(); ++i){
                if ((int)contigs[i].size() >= option.min_contig)
                    writer.Write(contigs[i], names[i]);
            }
            readTimehash_graph += MPI_Wtime();
            if(TIMING)printf("ReadSequence/writer time for kmer_size %i took %f seconds\n",kmer_size,readTimehash_graph);
            // scaffTime = -MPI_Wtime();
            // Scaffold(option.maxk, option.min_contig);
            // scaffTime += MPI_Wtime();
            // printf("Scaffold took %f seconds\n", scaffTime);
            string end_file = option.directory + "/end";
            fclose(OpenFile(end_file, "wb"));
            fflush(stdout);
        }  
        else{//these are all nodes that aren't the master node
            //save the file to parallel file system and send our rank to master to say we are done
            double writingTime = -MPI_Wtime();
            char filename [33];
            ofstream outfile;
            sprintf (filename,"%s/%d",option.directory.c_str(),kmer_size);
            outfile.open(filename, ios::binary | ios::out);
            outfile << hash_graph;
            outfile.close();
            writingTime += MPI_Wtime();
            if(TIMING)printf("Write time for rank %i and ksize %i took %f seconds\n",rank, kmer_size, writingTime);
            double sendTime = -MPI_Wtime();
            MPI_Send(&rank,1, MPI_INT,MASTER,0,MPI_COMM_WORLD);
            sendTime += MPI_Wtime();
            if(TIMING)printf("Sending from %i took %f seconds\n",rank,sendTime); 
        }
    }
    else{
        if(DEBUG)printf("%i Not being used\n",rank);
    }
    
    totalTime += MPI_Wtime();
    if(rank==MASTER)
        printf("Total Program time was %f seconds\n", totalTime);
    MPI_Finalize();
    return 0;
}
static void newBuildHashGraph(int kmer_size, HashGraph& hash_graph){
    char buffer [33];
    sprintf (buffer,"%d",kmer_size);
    BuildKmerFile(assembly_info, kmer_size, option.min_count, option.prefix_length, option.kmer_file()+buffer);
    ReadKmerFile(option.kmer_file()+buffer, hash_graph);    
    hash_graph.RefreshEdges();
    InsertInternalKmers(assembly_info, hash_graph, option.min_count);
}

static void UpdateHashGraph(int old_kmer_size,int kmer_size, HashGraph& hash_graph, int min_count){
    deque<Sequence> contigs;
    ReadSequence(option.contig_file(old_kmer_size), contigs);

    deque<Sequence> local_contigs;
    ReadSequence(option.local_contig_file(old_kmer_size), local_contigs);

    deque<Sequence> multi_contigs;
    ReadSequence(option.graph_file(old_kmer_size), multi_contigs);

    deque<Sequence> old_contigs;
    old_contigs.insert(old_contigs.end(), contigs.begin(), contigs.end());
    old_contigs.insert(old_contigs.end(), local_contigs.begin(), local_contigs.end());
    old_contigs.insert(old_contigs.end(), multi_contigs.begin(), multi_contigs.end());
    contigs.clear();
    local_contigs.clear();
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

    HashGraph tmp_hash_graph;
    tmp_hash_graph.swap(hash_graph);

    ContigGraph contig_graph(kmer_size, contigs, contig_infos);
    contigs.clear();
    contig_infos.clear();

    contig_graph.RemoveDeadEnd(option.min_contig);

    if (!option.is_no_bubble) {
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

    if (!option.is_no_coverage){
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

void LocalAssembly(int kmer_size, int new_kmer_size){
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

        cout << "edges " << scaffold_graph.num_edges(level) << endl;
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

void AddPairs(int level, ScaffoldGraph &scaffold_graph, const string &read_file, const string &align_file)
{
    ShortReadLibrary short_read_library;
    ReadLibrary(read_file, short_read_library);
    cout << "reads " << short_read_library.size() << endl;
    AlignReads(option.contig_file(option.maxk), short_read_library, align_file);

    EstimateDistance(align_file, median, sd);
    if (median < 0 || median != median || sd != sd || sd > 2*median){
        cout << "invalid insert distance" << endl;
        return;
    }

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
        if ((int)contigs[i].size() > option.min_contig){
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

