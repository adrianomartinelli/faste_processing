#include "extract_counts.h"

std::recursive_mutex m_fasta_file;
std::recursive_mutex m_counts;

ThreadStats::ThreadStats(StructFile &s) {
      const_mismatch = vector<unsigned long>(s.const_n_identifier.size());
      code_mismatch = vector<unsigned long>(s.code_n_identifier.size());

      const_match = vector<vector<unsigned long>>(s.const_n_identifier.size());
      for(int i = 0; i < const_match.size(); i++) {
        const_match[i].resize(s.const_n_identifier[i]);
      }

      code_match = vector<vector<unsigned long>>(s.code_n_identifier.size());
      for(int i = 0; i < code_match.size(); i++) {
        code_match[i].resize(s.code_n_identifier[i]);
      }
    };


//******Counts Class Function Definitions******//
size_t Counts::get_filesSize() { return file_size; }

vector<int> Counts::get_tulpel(unsigned long idx, StructFile &s, vector<int> loc, int n, unsigned long comb) {
  // idx%n_n, (idx%n_n)/n_n % n_n-1), ....
  vector<int> tulp;
  unsigned long coord;
  int k;

  if(idx > dim){
    cout<<"Index larger than dimension. Terminating"<<endl;
    exit(EXIT_FAILURE);
  }

  for (long i = n - 1; i >= 0; i--) {
    k = s.code_n_identifier[loc[i]];
    coord = idx % k;
    tulp.insert(tulp.begin(), coord);
    idx = idx / k;
  }
  return tulp;
}

void Counts::write_to_file(StructFile &s) {
  cout << "Writing selection files";
  string filename;
  vector<int> dim_tulpel;
  int stream_pre;
  int stream_post;

  for (unsigned long i = 0; i < comb_S; i++) {
    dim_tulpel = get_tulpel(i, s, S_idx, n_S, comb_S);

    filename = s.folder_name;
    filename.append("selection_");
    for_each(dim_tulpel.begin(), dim_tulpel.end(), [&filename](int x) {
      filename.append(to_string(x+1));
      filename.append("_");
    });
    filename.append(".txt");

    of.open(filename.c_str());
    if (!of.is_open()) {
      cout << "Cannot open file " << filename << " for writing." << endl;
      exit(EXIT_FAILURE);
    }
    if (of.is_open()) {
      of << "Count";
      for (int i = 0; i < n_B; i++) {
        of << "\tCode" << i+1;
      };
      of << endl;
      stream_pre = of.tellp();
      for (unsigned long j = 0; j < comb_B; j++) {
        // i*comb_B is the start index in the array for the ith selection, j loops over all possible bb combinations
        if (ignore_zero && C[i * comb_B + j]) {
          of << C[i * comb_B + j];
          print_vec(get_tulpel(i * comb_B + j, s, B_idx, n_B, comb_B), of);
          of << endl;
        } else if (!ignore_zero) {
          of << C[i * comb_B + j];
          print_vec(get_tulpel(i * comb_B + j, s, B_idx, n_B, comb_B), of);
          of << endl;
        }
      }

      stream_post = of.tellp();
      of.close();
    }

    if (stream_pre == stream_post && remove_empty_file) { // If nothing has been written to the file delete it.
      remove(filename.c_str());
    }
  }
  cout << ": \x1b[32mcompleted \x1b[0m" << endl;
}

unsigned long Counts::get_index(vector<unsigned int> &v) {
  unsigned long sum = 0;
  for (int i = 0; i < v.size(); i++) {
    sum += v[i] * cum_prod[i];
  }
  return sum;
  // assert(sum < dim);
};

void Counts::update_count(vector<unsigned int> &v) {
  unsigned long idx = get_index(v);
  guard g(m_counts);
  C[idx] += 1;
};

void Counts::print_vec(vector<int> v, ofstream &of) {
  for_each(v.begin(), v.end(), [this](int x) { this->of << "\t" << x+1; });
};

Counts::Counts(StructFile &s) {

  unsigned int tmp = 1;
  cum_prod.push_back(tmp);
  for (int i = s.code_n_identifier.size() - 1; i > 0; i--) {
    tmp *= s.code_n_identifier[i]; // M_N, M_N*M_{N-1}, ...
    cum_prod.insert(cum_prod.begin(), tmp);
  }

  dim = cum_prod[0] * s.code_n_identifier[0];
  C = new int[dim]();

  // Initialise B_idx, S_idx, n_S, n_B;
  for (int i = 0; i < s.is_selection.size(); i++) {
    if (s.is_selection[i] == 0) {
      B_idx.push_back(i);
      // n_B += s.code_n_identifier[i];
      n_B++;
      comb_B *= s.code_n_identifier[i];
    } else {
      S_idx.push_back(i);
      n_S++;
      comb_S *= s.code_n_identifier[i];
    }
  }

  // Get file size
  std::ifstream in(s.fasta_file.c_str(), std::ifstream::ate | std::ifstream::binary);
  file_size = in.tellg();
  in.close();
};

//******Function definitions******//
size_t fileposition(ifstream &ifs) {
  guard g(m_fasta_file);
  return ifs.tellg();
}

void thread_function(ifstream &inline_fasta, StructFile &s, ThreadStats &stats, Counts &counts, int i) {
  bool eof = false;
  string str_fasta;

  str_fasta = get_line(inline_fasta, eof);
  while (!eof) {
    analyse_sequence(str_fasta, s, stats, counts);
    stats.n_processed++;
    str_fasta = get_line(inline_fasta, eof);
  };
}

string get_line(ifstream &inline_fasta, bool &eof) {
  guard g(m_fasta_file);

  string str;
  if (!inline_fasta.eof()) {
    inline_fasta.ignore(1000, '\n'); // discar line starting with >
  } else {
    eof = true;
  }
  if (!inline_fasta.eof()) {
    getline(inline_fasta, str);
  } else {
    eof = true;
  }
  return str;
}

void analyse_sequence(string str_fasta, StructFile &s, ThreadStats &stats, Counts &counts) {
  bool const_mismatch = false;
  bool code_mismatch = false;
  bool match_flag = false;
  vector<unsigned int> idx;

  string str1;
  string str2;

  // Check const regions
  // Loop over all const lists
  for (int i = 0; i < s.const_sequences.size(); i++) {
    match_flag = false;

    for (int j = 0; j < s.const_sequences[i].size(); j++) {
      str1 = s.const_sequences[i][j];
      str2 = str_fasta.substr(s.const_location[i * 2], s.const_location[2 * i + 1]);

      if (compare_str(str1, str2)) {
        match_flag = true;
        stats.const_match[i][j]++; //increment match count for this const identifier
        break;
      }
    }
    if (match_flag) {
      match_flag = false;
    } else {
      const_mismatch = true;
      stats.const_mismatch[i]++;
      stats.n_discarded_seq++;
      break;
    }
  }

  // Check code regions
  if (!const_mismatch) {
    idx.clear();

    for (int i = 0; i < s.code_sequences.size(); i++) {
      match_flag = false;
      // Loop over all const sequences in the current list
      for (int j = 0; j < s.code_sequences[i].size(); j++) {
        str1 = s.code_sequences[i][j];
        str2 = str_fasta.substr(s.code_location[i * 2], s.code_location[2 * i + 1]);

        if (compare_str(str1, str2)) {
          match_flag = true;
          idx.push_back(j);
          stats.code_match[i][j]++; //increment match count for this const identifier
          break;
        }
      }
      if (match_flag) {
        match_flag = false;
      } else {
        code_mismatch = true;
        stats.code_mismatch[i]++;
        stats.n_discarded_seq++;
        break;
      }
    }
    if (!code_mismatch) {
      counts.update_count(idx);
      stats.n_analysed_seq++;
    }
  }
}

bool compare_str(string identifier, string substring) {
  return identifier == substring;
  // return !identifier.compare(substring);
}

void read_in(StructFile &s, ifstream &inline_struct) {

  int n; // sequence counter
  int start;
  int end;
  int code_type; // S: selection or B: building block
  string list_type;
  string list_file;
  string current_string;
  vector<string> current_sequences;
  ifstream inline_list;

  vector<vector<string>> *sequences;
  vector<int> *location;
  vector<int> *n_identifier;

  while (inline_struct >> start) {
    code_type = 0;
    // inline_struct >> start;
    inline_struct >> end;
    inline_struct >> list_type;
    inline_struct >> list_file;

    if (list_type == "C") {
      sequences = &s.const_sequences;
      location = &s.const_location;
      n_identifier = &s.const_n_identifier;
    } else if (list_type == "S" || list_type == "B") {
      sequences = &s.code_sequences;
      location = &s.code_location;
      n_identifier = &s.code_n_identifier;
      if (list_type == "S") {
        code_type = 1;
      }
      s.is_selection.push_back(code_type);
    } else {
      std::cerr << "Error: Invalid list specification. Use either S for selections code lists, B for building block code lists or C for const regions" << endl;
      exit(EXIT_FAILURE);
    }

    location->push_back(start-1);
    location->push_back(end - start + 1);

    inline_list.open(list_file.c_str());
    if (!inline_list.is_open()) {
      cout << "Cannot open file " << list_file << " for reading." << endl;
      exit(EXIT_FAILURE);
    }
    if (inline_list.is_open()) {
      cout << "  Reading in " << list_file << endl;
      current_sequences.clear();
      n = 0;

      while (inline_list >> current_string) {
        current_sequences.push_back(current_string);
        n++;
      }
      sequences->push_back(current_sequences);
      n_identifier->push_back(n);
      inline_list.close();
    }
  }
}

void summary(StructFile &s, Statistics &stats, ostream &stream) {
  stream << "---------- Overview ----------" << endl;
  stream << "Analysed file: " << s.fasta_file << endl;
  stream << "Total number of sequences received:\t" << stats.n_processed << endl;
  stream << "Total number of sequences processed:\t" << stats.n_analysed_seq + stats.n_discarded_seq << " (" << stats.percent_processed << "%)" << endl;
  stream << "Total number of sequences analysed:\t" << stats.n_analysed_seq << " (" << stats.percent_analysed << " %)" << endl;
  stream << "Total number of sequences discarded:\t" << stats.n_discarded_seq << " (" << stats.percent_discarded << " %)" << endl;
  stream << endl;
  stream << "---------- Mismatch Report ----------" << endl;
  stream << "Constant Region";
  string buff = "\t";
  unsigned long sum = 0;
  for (int i = 0; i < stats.const_mismatch.size(); i++) {
    stream << "\t" << i;
    buff += "\t" + to_string(stats.const_mismatch[i]);
    sum += stats.const_mismatch[i];
  }
  stream << endl;
  stream << buff << endl << "Total: " << sum << endl << endl;

  stream << "Code Region";
  buff = "\t";
  sum = 0;
  for (int i = 0; i < stats.code_mismatch.size(); i++) {
    stream << "\t" << i;
    buff += "\t" + to_string(stats.code_mismatch[i]);
    sum += stats.code_mismatch[i];
  }
  stream << endl;
  stream << buff << endl << "Total: " << sum << endl << endl;

  stream << "---------- Struct File ----------" << endl;
  for (int i = 0; i < s.const_sequences.size(); i++) {
    stream << "Constant Region " << i << endl
           << "Number of identifiers: " << s.const_n_identifier[i] << endl
           << "Start: " << s.const_location[2 * i] << "\tEnd: " << s.const_location[2 * i] + s.const_location[2 * i + 1] - 1 << "\tLength: " << s.const_location[2 * i + 1] << "\tType: C" << endl;
    //printn_vec(s.const_sequences[i], stream);
    sum = 0;
    for(int j = 0; j < s.const_sequences[i].size(); j++){
      sum += stats.const_match[i][j];
      stream << j+1 << '\t' << s.const_sequences[i][j] << "\t" << stats.const_match[i][j] <<endl;
      }
    stream << "Total\t" << sum << endl;
    stream << endl;
  }

  for (int i = 0; i < s.code_sequences.size(); i++) {
    string type = "B";
    if (s.is_selection[i]) {
      type = "S";
    }
    stream << "Code Region " << i << endl
           << "Number of identifiers: " << s.code_n_identifier[i] << endl
           << "Start: " << s.code_location[2 * i] << "\tEnd: " << s.code_location[2 * i] + s.code_location[2 * i + 1] + 1 << "\tLength: " << s.code_location[2 * i + 1] << "\tType: " << type << endl;
    //printn_vec(s.code_sequences[i], stream);
    sum = 0;
    for(int j = 0; j < s.code_sequences[i].size(); j++){
      sum += stats.code_match[i][j];
      stream << j+1 << '\t' << s.code_sequences[i][j] << "\t" << stats.code_match[i][j] <<endl;
    }
    stream << "Total\t" << sum << endl;

    stream << endl;
  }
}

template <class T> void print_vec(vector<T> v) {
  for_each(v.begin(), v.end(), [](T x) { cout << x << " "; });
}