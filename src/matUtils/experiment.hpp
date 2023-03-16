#include "convert.hpp"
#include "select.hpp"
#include <string>
#include <random>
#include <regex>
#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>

struct ances_sample_list {
   Mutation_Annotated_Tree::Node * ancestor_node; 
   std::vector<Mutation_Annotated_Tree::Node*> * sample_nodes;
};

struct sample_ances_list {
   Mutation_Annotated_Tree::Node * sample; 
   std::vector<Mutation_Annotated_Tree::Node*> * ances;
};

struct sample_read_pair {
   Mutation_Annotated_Tree::Node * sample;
   std::string read;
};

struct pos_misread {
   int pos;
   std::string read;
   bool used;
};

struct read_mut_pair {
   std::string read_name;
   char mut_nuc;
};

struct read_info {
   std::string read;
   std::vector<MAT::Mutation> mutations;
   int start;
   int end;
};

struct parsimony {
   std::vector<MAT::Mutation> p_node_par;
   MAT::Node* curr_node; 
};

struct min_parsimony {
   std::vector<size_t> idx_list;
   std::vector<std::vector<MAT::Mutation>> par_list;
   std::vector<bool> is_sibling_list;
};

po::variables_map parse_extract_command(po::parsed_options parsed);
void simulate_and_place_reads(po::parsed_options parsed);

void read_vcf(std::string, std::vector<struct read_info*> &);
void place_reads(const std::vector<MAT::Node*> &, const std::vector<struct read_info*> &, tbb::concurrent_hash_map<size_t, struct min_parsimony> &, std::string);
