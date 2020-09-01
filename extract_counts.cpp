/*READ ME
Requirements:
Not new line at the end of struct file
Selections bar codes before building block bar codes

Run the following code snipped in your terminal to compile a executable file
g++ -O3 -std=c++11 -o extract_counts.o extract_counts.cpp extract_counts_fun.cpp


Run executable with:
./extract_counts_final.o

Rename file in folder with:
for j in *_.txt; do mv $j 'selection_'$j; done
*/

#include <algorithm> //for_each
#include <chrono>    //Wall time
#include <ctime>     //CPU time
#include <functional> //std::function
#include <sys/stat.h>
#include <thread> //Multithreading
#include "extract_counts.h"

using namespace std;

int main(int argc, char *const argv[]) {

  // Measure wall time and and CPU time
  chrono::time_point<chrono::system_clock> start, end, start_read, end_read, start_fasta, end_fasta, start_write, end_write;
  clock_t start_CPU, end_CPU, start_read_CPU, end_read_CPU, start_fast_CPU, end_fasta_CPU, start_write_CPU, end_write_CPU;
  start = chrono::system_clock::now();
  start_CPU = clock();

  // Files from which we read in
  // Files from which we read in
  string structure_file;
  if(argc == 1){
    cout << "Specify structure file: ";
    cin >> structure_file;
    }
  else{
    structure_file = argv[1];
  }
  bool summary_to_file = true; // write summary file;

  // Stores all information from struct file
  StructFile s;

  // Multithreading maximum available number cores
  int cores = thread::hardware_concurrency();

  ifstream inline_struct(structure_file.c_str());
  if (!inline_struct) {
    cout << "Cannot open file " << structure_file << " for reading." << endl;
    exit(EXIT_FAILURE);
  }
  if (inline_struct.is_open()) {
    start_read = chrono::system_clock::now();
    start_read_CPU = clock();

    cout << "Reading in structure file: " << structure_file << endl;
    inline_struct >> s.fasta_file;
    inline_struct >> s.folder_name;

    // Read in sequence lists
    read_in(s, inline_struct);

    inline_struct.close();
    cout << "Reading in struct file: \x1b[32mcompleted \x1b[0m" << endl;

    end_read = chrono::system_clock::now();
    end_read_CPU = clock();
  }

  // Creat folder to store selection files
  const int dir_err = mkdir(s.folder_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  // Data data(s.code_n_identifier);
  Counts counts(s);

  // Statistic objects for each thread
  vector<ThreadStats> v_threadStats(cores, {s});

  // Process fasta file
  ifstream infile_fasta_file(s.fasta_file.c_str());
  if (!infile_fasta_file) {
    cout << "Cannot open file " << s.fasta_file << " for reading." << endl;
    exit(EXIT_FAILURE);
  }
  if (infile_fasta_file.is_open()) {
    start_fasta = chrono::system_clock::now();
    start_fast_CPU = clock();

    cout << "Processing fasta file" << endl;
    vector<thread> v_threads;
    for (int i = 0; i < cores; i++) {
      v_threadStats[i].thread_id = i;
      v_threads.push_back(std::thread(thread_function, ref(infile_fasta_file), ref(s), ref(v_threadStats[i]), ref(counts), i));
    }

    ProgressBar<size_t> pb(0, counts.get_filesSize());
    bool progress = true;
    size_t tell;
    while (progress) {
      tell = fileposition(infile_fasta_file);
      progress = pb.progress(tell);
      this_thread::sleep_for(chrono::milliseconds(1000));
    }
    cout << endl;

    for_each(v_threads.begin(), v_threads.end(), [](thread &t) { t.join(); });
    infile_fasta_file.close();
    cout << "Processing fasta file: \x1b[32mcompleted \x1b[0m" << endl;

    end_fasta = chrono::system_clock::now();
    end_fasta_CPU = clock();
  }

  // Initialise Stats objects
  Statistics stats(v_threadStats);

  // Write counts to files
  start_write = chrono::system_clock::now();
  start_write_CPU = clock();
  counts.write_to_file(s);
  end_write = chrono::system_clock::now();
  end_write_CPU = clock();

  if (summary_to_file) {
    cout << "Writing summary file: ";
    ofstream stream;

    string filename = s.folder_name;
    filename.append("summary.txt");

    stream.open(filename.c_str());
    if (!stream.is_open()) {
      stream << "Cannot open file " << filename << " for writing." << endl;
      exit(EXIT_FAILURE);
    }
    if (stream.is_open()) {
      summary(s, stats, stream);
    }
    stream.close();
    cout << " \x1b[32mcompleted \x1b[0m" << endl << endl;
  }
  summary(s, stats, cout); // Output summary to console

  end = chrono::system_clock::now();
  end_CPU = clock();

  double duration_CPU, duration_read_CPU, duration_fasta_CPU, duration_write_CPU;
  chrono::duration<double> duration, duration_read, duration_fasta, duration_write;

  duration_CPU = 1.0 * (end_CPU - start_CPU) / CLOCKS_PER_SEC;
  duration_read_CPU = 1.0 * (end_read_CPU - start_read_CPU) / CLOCKS_PER_SEC;
  duration_fasta_CPU = 1.0 * (end_fasta_CPU - start_fast_CPU) / CLOCKS_PER_SEC;
  duration_write_CPU = 1.0 * (end_write_CPU - start_write_CPU) / CLOCKS_PER_SEC;
  duration = (end - start);
  duration_read = (end_read - start_read);
  duration_fasta = (end_fasta - start_fasta);
  duration_write = (end_write - start_write);

  cout << "---------- Computation Time ----------" << endl;
  cout << "Elapsed Time: " << duration.count() << endl;
  cout << "  Read Files: " << duration_read.count() << endl;
  cout << "  Process Sequences: " << duration_fasta.count() << endl;
  cout << "  Write Files: " << duration_write.count() << endl;
  cout << "CPU Time: " << duration_CPU << endl;
  cout << "  Read Files: " << duration_read_CPU << endl;
  cout << "  Process Sequences: " << duration_fasta_CPU << endl;
  cout << "  Write Files: " << duration_write_CPU << endl;

  return 0;
}
