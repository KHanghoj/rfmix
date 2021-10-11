/* RFMIX v2.XX - Local Ancestry and Admixture Analysis
   Bustamante Lab - Stanford School of Medicine
   (c) 2016 Mark Hamilton Wright

   This program is licensed for academic research use only
   unless otherwise stated. Contact cdbadmin@stanford.edu for
   commercial licensing options.

   Academic and research users should cite Brian Maples'
   paper describing RFMIX in any publication using RFMIX
   results. Citation is printed when the program is started. */

#include <fstream>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <unistd.h>
#include <stdint.h>
#include <limits.h>

#include "kmacros.h"
#include "cmdline-utils.h"
#include "vcf.h"
#include "s-sample.h"
#include "s-subpop.h"

typedef struct {
  char *vcf_fname;
  char *sample_map_fname;
  char *genetic_fname;
  char *output_basename;
  char *chromosome;
  char *mating_fname;
  int n_generations;
  int ril;
  double growth_rate;
  char *random_seed_str;
  int32_t random_seed;
  double phase_switch_rate;
} opts_t;

opts_t opts;

static option_t options[] = {
  { 'f', "vcf", &opts.vcf_fname, OPT_STR, 1, 1,
    "Name of input TRALA VCF file" },
  { 'x', "matingfile", &opts.mating_fname, OPT_STR, 1, 1,
    "Name of file from sampler.py" },
  { 'm', "sample-map", &opts.sample_map_fname, OPT_STR, 1, 1,
    "Sample subpop mapping file - also selects which samples will be used for simulation" },
  { 'g', "genetic-map", &opts.genetic_fname, OPT_STR, 1, 1,
    "Genetic map file (required)" },
  { 'o', "output-basename", &opts.output_basename, OPT_STR, 1, 1,
    "Basename (prefix) for output files (required)" },
  { 0, "growth-rate", &opts.growth_rate, OPT_DBL, 0, 1,
    "Growth rate of population per generation (1 = no growth, 1.5 = 50% increase per generation, etc.)" },  
  { 'c', "chromosome", &opts.chromosome, OPT_STR, 1, 1,
    "Chromosome to select from the VCF file" },
  {  0 , "make-rils", &opts.ril, OPT_FLAG, 0, 0,
     "After first generation of random mating, make recombinant inbred lines by selfing" },
  { 'p', "phasing-switch", &opts.phase_switch_rate, OPT_DBL, 0, 1,
    "introduce phasing switches at this rate. verification data is not switched" },
  
  { 0, "random-seed", &opts.random_seed_str, OPT_STR, 0, 1,
    "Seed value for random number generation - integer value (maybe specified in"
    "hexadecimal by preceeding with 0x), or the string \"clock\" to seed with "
    "the current system time." },
  { 0, NULL, NULL, 0, 0, 0, NULL }
};

static void init_options(void) {
  opts.vcf_fname = NULL;
  opts.sample_map_fname = NULL;
  opts.mating_fname = NULL;
  opts.genetic_fname = NULL;
  opts.output_basename = NULL;
  opts.growth_rate = 1.0;
  opts.chromosome = NULL;
  opts.ril = 0;
  opts.random_seed_str = (char *) "0xDEADBEEF";
  opts.phase_switch_rate = 0.0;
}

static void verify_options(void) {

  if (opts.vcf_fname == NULL) {
    fprintf(stderr,"\nSpecify VCF input file for source data with -f option\n\n");
    exit(-1);
  }
  if (opts.mating_fname == NULL) {
    fprintf(stderr,"\nSpecify matingfile -x option\n\n");
    exit(-1);
  }
  if (opts.output_basename == NULL) {
    fprintf(stderr,"\nSpecify the output basename (prefix) with -o option\n\n");
    exit(-1);
  }
  if (opts.genetic_fname == NULL) {
    fprintf(stderr,"\nA genetic map is required, specify with -g option\n\n");
    exit(-1);
  }
  if (opts.sample_map_fname == NULL) {
    fprintf(stderr,"\nSpecify sample to subpopulation mapping file with -m option\n\n");
    exit(-1);
  }

  if (opts.chromosome == NULL) {
    fprintf(stderr,"\nSpecify chromosome to select from VCF input with -c option\n\n");
    exit(-1);
  }
  if (opts.phase_switch_rate < 0. || opts.phase_switch_rate > 1.0) {
    fprintf(stderr,"\nPhasing switch rate (-p) must be between 0 and 1\n\n");
    exit(-1);
  }
  
  if (strcmp(opts.random_seed_str, "clock") == 0) {
    opts.random_seed = time(NULL);
  } else {
    opts.random_seed = strtod(opts.random_seed_str,0);
  }
  srand(opts.random_seed);  
}

#if 0
static int load_subpop_map(char **subpops, S_Sample **samples, char *fname) {
  f = fopen(fname, "r");
  if (f == NULL) {
    fprintf(stderr,"Can't open input file %s (%s)\n", fname, strerror(errno));
    exit(-1);
  }

  while(fgets(buf, 8192, f) != NULL) {
    CHOMP(buf);
    if (buf[0] == 0 || buf[0] == '#') continue;

    char *p = buf;
    char *sample_id = strsep(&p, "\t");
    char *subpop_name = strsep(&p, "\t");

    int subpop_idx = 0;
    while(subpop_idx < n_subpops && strcmp(subpop_map[subpop_idx], subpop_name) != 0) subpop_idx++;
    if (subpop_idx == n_subpops) {
      if (n_subpops % SUBPOP_ALLOC_STEP == 0)
	RA(subpop_map, sizeof(char *)*(n_subpops + SUBPOP_ALLOC_STEP), char *);
      subpop_map[n_subpops] = strdup(subpop_name);
      n_subpops++;
    }

    for(int i=0; i < n_samples; i++) {
    }
  }
}
#endif

#define BUF_SIZE (8192)
static void load_sample_subpop_map(char *fname) {
  char buf[BUF_SIZE];
  
  FILE *f = fopen(fname, "r");
  if (f == NULL) {
    fprintf(stderr,"Can't open input file %s (%s)\n", fname, strerror(errno));
    exit(-1);
  }

  while(fgets(buf, 8192, f) != NULL) {
    CHOMP(buf);
    if (buf[0] == 0 || buf[0] == '#') continue;

    char *p = buf;
    char *sample_name = strsep(&p, "\t");
    char *subpop_name = strsep(&p, "\t");
    
    Subpop *s = Subpop::lookup_subpop(subpop_name);
    if (s == NULL) s = new Subpop(subpop_name);

    s->add_sample(sample_name);
  }
  
  fclose(f);
  Subpop::set_ordering();
}

std::vector<std::vector<int>> read_mating_file(char *fname){
  
  std::vector<std::vector<int>> res;
  std::vector<int> temp;
  std::string line, word;
  std::istringstream ss;
  
  std::ifstream fh(fname);
  if(!fh.good()){
    fprintf(stderr,"Can't open input file %s\n", fname);
    exit(-1);
  }

  while(std::getline(fh, line)) {
    // fprintf(stderr, "%s\n", line.c_str());
    ss.clear();
    ss.str(line);
    while(ss >> word){
      temp.push_back(std::stoi(word));
    }
    res.push_back(temp);
    temp.clear();
  }
  fh.close();
  
  return(res);
}
 

int main(int argc, char *argv[]) {

  init_options();
  cmdline_getoptions(options, argc, argv);
  verify_options();


  std::vector<std::vector<int>> mating_pattern = read_mating_file(opts.mating_fname);
  int n_generations = mating_pattern.size();
#if 0
  fprintf(stderr, "%ld\n",mating_pattern.size());
  // remember that n_ind == size()//2;
  for (size_t i=0; i<mating_pattern.size();i++){
    for (size_t j=0; j<mating_pattern[i].size();j++){
      fprintf(stderr, "%d", mating_pattern[i][j]);
    }
    fprintf(stderr, "\n");
  }
  exit(0);
#endif
  
  fprintf(stderr,"\nLoading VCF %s chromosome %s ... ", opts.vcf_fname, opts.chromosome);
  VCF *vcf = new VCF(opts.vcf_fname, opts.chromosome);
  load_sample_subpop_map(opts.sample_map_fname);

  GeneticMap *genetic_map = new GeneticMap();
  genetic_map->load_map(opts.genetic_fname, opts.chromosome);
  
  vcf->load_snps(opts.chromosome, genetic_map);
  vcf->load_haplotypes(opts.chromosome);

  int n_samples = 0;
  S_Sample **parents = new S_Sample*[vcf->n_samples];
  for(int i=0; i < vcf->n_samples; i++) {
    Subpop *s = Subpop::lookup_sample_subpop(vcf->samples[i].sample_id);
    if (s != NULL)
      parents[n_samples++] = new S_Sample(vcf->samples[i].sample_id, s->idx, vcf->snps, vcf->n_snps,
					vcf->samples[i].haplotype[0], vcf->samples[i].haplotype[1]);
  }
  fprintf(stderr,"\nChromosome %s loaded - %d SNPs, %d samples\n", opts.chromosome, vcf->n_snps, n_samples);
  vcf->DropHaplotypes();
  
  int last_size = n_samples;
  for(int g=0; g < n_generations; g++) {
    int next_size = mating_pattern[g].size() / 2;
    if (isatty(2)) {
      fprintf(stderr,"\rBreeding generation %4d/%4d (%1.1f%%) - %6d diploid individuals     ",
	      g+1, n_generations, (g+1)/(double) n_generations * 100., next_size);
    }

    // int i;
    // int *s = new int[last_size];
    // for(i=0; i < last_size; i++)
    //   s[i] = i;

    // fprintf(stderr, "%d %d %d\n", g, last_size, next_size);        
    S_Sample **children = new S_Sample*[next_size];
    int mating_idx = 0;
    for (int child_idx=0; child_idx<next_size; child_idx++, mating_idx+=2){
      int p1_idx = mating_pattern[g][mating_idx];
      int p2_idx = mating_pattern[g][mating_idx+1];
      // fprintf(stderr, "%d %d %d %d %d\n", child_idx, mating_idx, mating_idx+1, p1_idx, p2_idx);
      children[child_idx] = new S_Sample(parents[p1_idx], parents[p2_idx]);
    }
    // for(int i=0; i < next_size; i++) {
    //   if (i % last_size == 0) {
    //     for(int k=0; k < last_size; k++) {
    //       int j = rand()/(RAND_MAX + 1.0) * last_size;
    //       int tmp = s[k];
    //       s[k] = s[j];
    //       s[j] = tmp;
    //     }
    //   }
      
      // int p1_idx, p2_idx;
      // if (opts.ril == 1 && g > 0) {
      //   p1_idx = p2_idx = s[i % last_size];
      // } else {
      //   p1_idx = s[i % last_size];
      //   p2_idx = s[(i+1) % last_size];
      // }
      
      // children[i] = new S_Sample(parents[p1_idx], parents[p2_idx]);
    // }
    // delete[] s;
    
    for(int i=0; i < last_size; i++)
      delete parents[i];
    delete[] parents;

    last_size = next_size;
    parents = children;
  }
  fprintf(stderr,"done.\nFinal simulated population size is %d diploid individuals\n", last_size);

  fprintf(stderr,"Writing output... ");
  char vcf_out_fname[strlen(opts.output_basename) + strlen(".query.vcf") + 1];
  sprintf(vcf_out_fname,"%s.query.vcf", opts.output_basename);
  char result_fname[strlen(opts.output_basename) + strlen(".result") + 1];
  sprintf(result_fname,"%s.result", opts.output_basename);
  
  FILE *vf = fopen(vcf_out_fname, "w");
  if (vf == NULL) {
    fprintf(stderr,"Can't open VCF output file %s (%s)\n", vcf_out_fname, strerror(errno));
    exit(-1);
  }
  FILE *rf = fopen(result_fname, "w");
  if (rf == NULL) {
    fprintf(stderr,"Can't open result output file %s (%s)\n", result_fname, strerror(errno));
    exit(-1);
  }

  int phase_switched[last_size];
  fprintf(vf,"##fileformat=VCFv4.1\n");
  fprintf(vf,"##source=%s (rfmix v2)\n", argv[0]);
  fprintf(vf,"##FORMAT=<ID=GT,Number=1,Type=String,Descripton=\"Phased Genotype\">\n");
  fprintf(vf,"##contig=<ID=%s>\n", opts.chromosome);
  fprintf(vf,"#CHROM\tPOS\tID\tREF\tVAR\tQUAL\tFILTER\tINFO\tFORMAT");
  fprintf(rf,"chm\tpos");
  for(int j=0; j < last_size; j++) {
    phase_switched[j] = 0;
    fprintf(vf,"\tsim%d", j);
    fprintf(rf,"\tsim%d.0\tsim%d.1", j, j);
    // fprintf(vf,"\t%s", parents[j]->sample_id);
    // fprintf(rf,"\t%s.0\t%s.1", parents[j]->sample_id, parents[j]->sample_id);
  }
  fprintf(vf,"\n");
  fprintf(rf,"\n");

  int i;
  for(i=0; i < vcf->n_snps; i++) {
    if (((i & 0xFFF) == 0) && isatty(2)) {
      fprintf(stderr,"\rWriting output...  %11d/%11d (%1.1f%%) SNPs     ", i, vcf->n_snps,
	      i/(double) vcf->n_snps * 100.);
    }
    fprintf(vf,"%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT", opts.chromosome, vcf->snps[i].pos,
	    vcf->snps[i].snp_id, vcf->snps[i].ref, vcf->snps[i].alt);
    fprintf(rf,"%s\t%d", opts.chromosome, vcf->snps[i].pos);

    for(int j=0; j < last_size; j++) {
      if (rand()/(RAND_MAX + 1.0) < opts.phase_switch_rate) {
	phase_switched[j] = phase_switched[j] ? 0 : 1;
      }
    }
    for(int j=0; j < last_size; j++) {
      fputc('\t',vf);
      if (phase_switched[j]) {
	fprintf(vf,"%c|%c",parents[j]->haplotype[1][i] + '0',
		parents[j]->haplotype[0][i] + '0');
      } else {
	fprintf(vf,"%c|%c",parents[j]->haplotype[0][i] + '0',
		parents[j]->haplotype[1][i] + '0');
      }
      //      fprintf(vf,"\t%d|%d", parents[j]->haplotype[0][i], parents[j]->haplotype[1][i]);
      fputc('\t',rf);
      fputc(parents[j]->subpop[0][i] + '1', rf);
      fputc('\t',rf);
      fputc(parents[j]->subpop[1][i] + '1', rf);
      //      fprintf(rf,"\t%d\t%d", parents[j]->subpop[0][i] + 1, parents[j]->subpop[1][i] + 1);
    }
    fprintf(vf,"\n");
    fprintf(rf,"\n");
  }
  if (isatty(2))
    fprintf(stderr,"\r");
  else
    fprintf(stderr,"\n");
  fprintf(stderr,"Writing output...  %11d/%11d (%1.1f) SNPs     done.\n", i, vcf->n_snps,
	  i/(double) vcf->n_snps * 100.);
  
  fclose(vf);
  fclose(rf);

  delete[] parents;
  delete vcf;
  delete genetic_map;
  return 0;
}


