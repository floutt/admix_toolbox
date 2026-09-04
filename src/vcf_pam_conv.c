#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <dotgeno.h>

#define PAM_STR_BUF_EXTRA 6
#define VCF_SAMP_ADJ(x) (x - 9)

char** str_split(char* str, char delim, size_t* n_elems) {
	size_t str_len = strlen(str);
	size_t nc_buf = str_len > 0;
	// get num chars
	for(size_t i = 0; i < str_len; i++) {
		if(str[i] == delim) {
			nc_buf += 1;
		}
	}
	size_t* col_lens = (size_t*)malloc(nc_buf * sizeof(size_t)); 
	size_t cur_len = 0;
	size_t cur_col = 0;
	for(size_t i = 0; i < str_len; i++) {
		if(!(str[i] == delim)) {
			cur_len += 1;
		} else {
			col_lens[cur_col] = cur_len;
			cur_len = 0;
			cur_col += 1;
		}
	}
	col_lens[cur_col] = cur_len;
	char** out_arr = (char**)malloc(nc_buf * sizeof(char*));
	for(size_t i = 0; i < nc_buf; i++) {
		out_arr[i] = (char*)malloc((col_lens[i] + 1) * sizeof(char));
		out_arr[i][col_lens[i]] = '\0';
	}
	size_t subber = 0;
	cur_col = 0;
	for(size_t i = 0; i < str_len; i++) {
		if((str[i] == delim)) {
			cur_col += 1;
			subber = i + 1;
		} else {
			out_arr[cur_col][i - subber] = str[i];
		}
	}
	free(col_lens);
	*n_elems = nc_buf;
	return out_arr;
}

void write_snp_and_ind(FILE* f_vcf, char* snp_file, char* ind_file, char** ind_pop, char** ind_sex, size_t ind_len) {
	// go to beginning of file
	fseek(f_vcf, 0, SEEK_SET);
	
	// write ind
	char chr[7] = "#CHROM";
	char* line = NULL;
	size_t size = 0;
	ssize_t nread;
	while ((nread = getline(&line, &size, f_vcf)) != -1) {
		if(line[0] != '#') {
			fprintf(stderr, "ERROR: could not find #CHROM header in VCF!");
			exit(EXIT_FAILURE);
		} else if(strncmp(line, chr, 6) == 0) {
			break;
		}
	}
	line[nread - 1] = '\0';
	size_t n_elems;
	char** elems = str_split(line, '\t', &n_elems);
	if(VCF_SAMP_ADJ(n_elems) != ind_len) {
		fprintf(stderr, "ERROR: invalid number of individuals provided. Expected %zu, got %zu\n", VCF_SAMP_ADJ(n_elems), ind_len);
		exit(EXIT_FAILURE);
	}

	FILE* f_ind = fopen(ind_file, "w+");
	for(size_t i = 9; i < n_elems; i++) {
		char* sex = "U";
		if(ind_sex) { sex = ind_sex[VCF_SAMP_ADJ(i)]; }
		fprintf(f_ind, "%s\t%s\t%s\n", elems[i], sex, ind_pop[VCF_SAMP_ADJ(i)]);
	}
	fclose(f_ind);
	for(size_t i = 0; i < n_elems; i++) {
		free(elems[i]);
	}
	free(elems);
	// write snp
	FILE* f_snp = fopen(snp_file, "w+");
	while ((nread = getline(&line, &size, f_vcf)) != -1) {
		line[nread - 1] = '\0';
		elems = str_split(line, '\t', &n_elems);
		char* chrom = elems[0];
		char* pos = elems[1];
		char* snp_id = elems[2];
		char* ref = elems[3];
		char* alt = elems[4];
		// skip if multiallelic
		bool is_multi_allele = false;
		size_t len_alt = strlen(alt);
		for(size_t i = 0; i < len_alt; i++) {
			if(alt[i] == ',') {
				is_multi_allele = true;
				break;
			}
		}
		if(is_multi_allele) {
			for(size_t i = 0; i < n_elems; i++) {
				free(elems[i]);
			}
			free(elems);
			continue;
		}

		// if not rsid then change snp_id
		if(strncmp("rs", snp_id, 2) != 0) {
			snp_id = (char*)malloc(sizeof(char) * (strlen("SNP") + strlen(chrom) + strlen(pos) + 3));
			sprintf(snp_id, "snp_%s_%s", chrom, pos);
		}
		fprintf(f_snp, "%s\t%s\t%s\t%s\t%s\t%s\n", snp_id, chrom, "0.0", pos, ref, alt);
		if(strncmp("snp", snp_id, 3) == 0) { free(snp_id); } 
		for(size_t i = 0; i < n_elems; i++) {
			free(elems[i]);
		}
		free(elems);
	}
	fclose(f_snp);
	// free
	free(line);
}

void write_geno(FILE* f_vcf, snp_data* snp, ind_data* ind, char* geno_file) {
	// go to beginning of file
	fseek(f_vcf, 0, SEEK_SET);
	
	char* line = NULL;
	size_t size = 0;
	ssize_t nread;

	pam_file_writer pfw = pam_file_writer_init(geno_file, snp, ind);
	write_pam_header(&pfw, snp, ind);
	while ((nread = getline(&line, &size, f_vcf)) != -1) {
		line[nread - 1] = '\0';
		if(line[0] == '#') { continue; }
		size_t n_elems;
		char** elems = str_split(line, '\t', &n_elems);
		if(VCF_SAMP_ADJ(n_elems) != ind->length) {
			fprintf(stderr, "ERROR: invalid number of individuals provided in VCF reecord. Expected %zu, got %zu\n", VCF_SAMP_ADJ(n_elems), ind->length);
			exit(EXIT_FAILURE);
		}

		bool is_multi_allele = false;
		char* alt =  elems[4];
		size_t len_alt = strlen(alt);
		for(size_t i = 0; i < len_alt; i++) {
			if(alt[i] == ',') {
				is_multi_allele = true;
				break;
			}
		}
		if(is_multi_allele) {
			for(size_t i = 0; i < n_elems; i++) {
				free(elems[i]);
			}
			free(elems);
			continue;
		}
		// get GT num
		size_t n_fmt;
		size_t gt_pos;
		bool found_gt = false;
		char** fmt_str_split = str_split(elems[8], ':', &n_fmt);
		for(size_t i = 0; i < n_fmt; i++) {
			if(strcmp(fmt_str_split[i], "GT") == 0) {
				found_gt = true;
				gt_pos = i;
				break;
			}
		}
		for(size_t i = 0; i < n_fmt; i++) {
			free(fmt_str_split[i]);
		}
		free(fmt_str_split);

		if(!found_gt) {
			fprintf(stderr, "ERROR: GT not found in row.\n");
			exit(EXIT_FAILURE);
		}

		uint8_t* record = (uint8_t*)malloc(sizeof(uint8_t) * ind->length);
		for(size_t i = 9; i < n_elems; i++) {
			char* samp_str = elems[i];
			char** samp_str_split = str_split(samp_str, ':', &n_fmt);
			char* gt_str = samp_str_split[gt_pos];
			uint8_t dosage;
			if((strcmp("0/0", gt_str) == 0) || (strcmp("0|0", gt_str) == 0)) {
				dosage = 2;
			} else if((strcmp("1/1", gt_str) == 0) || (strcmp("1|1", gt_str) == 0)) {
				dosage = 0;
			} else if((strcmp("1/0", gt_str) == 0) || (strcmp("1|0", gt_str) == 0) || (strcmp("0/1", gt_str) == 0) || (strcmp("0|1", gt_str) == 0)) {
				dosage = 1;
			} else if((strcmp("./.", gt_str) == 0) || (strcmp(".", gt_str) == 0)) {
				dosage = NAN_VAL;
			} else {
				fprintf(stderr, "ERROR: invalid value '%s' in sample column of VCF.\n", gt_str);
				exit(EXIT_FAILURE);
			}
			record[VCF_SAMP_ADJ(i)] = dosage;

			// free samp_str_split
			for(size_t j = 0; j < n_fmt; j++) {
				free(samp_str_split[j]);
			}
			free(samp_str_split);
		}
		write_pam_record(&pfw, record);
		free(record);
		for(size_t i = 0; i < n_elems; i++) {
			free(elems[i]);
		}
		free(elems);
	}
	free(line);
	close_pam_file_writer(&pfw);
}

void write_vcf(pam_file_reader* pfr, char* out_vcf, snp_data* snp, ind_data* ind) {
	FILE* f_out = fopen(out_vcf, "w+");
	uint8_t* record;
	fprintf(f_out, "##fileformat=VCFv4.2\n");
	fprintf(f_out, "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
	fprintf(f_out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	char* last_chr = NULL;
	// assumes SNP is sorted by chr
	for(size_t i = 0; i < snp->length; i++) {
		if((last_chr == NULL) || (strcmp(last_chr, snp->chr[i]) != 0)) {
			fprintf(f_out, "##contig=<ID=%s>\n", snp->chr[i]);
			last_chr = snp->chr[i];
		}
	}
	fprintf(f_out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
	for(size_t i = 0; i < ind->length; i++) {
		fprintf(f_out, "%s", ind->ind_id[i]);
		if(i == (ind->length - 1)) { fprintf(f_out, "\n"); }
		else { fprintf(f_out, "\t"); }
	}

	while(record = read_pam_record(pfr)) {
		size_t idx = pfr->idx - 1;
		fprintf(f_out, "%s\t%" PRId64 "\t%s\t%s\t%s\t%s\t%s\t%s", snp->chr[idx], snp->pos[idx], snp->var_id[idx], snp->ref[idx], snp->alt[idx], "999", "PASS", ".", "GT");
		for(size_t i = 0; i < pfr->n_ind; i++) {
			switch(record[i]) {
				case 0:
					fprintf(f_out, "\t1/1");
					break;
				case 1:
					fprintf(f_out, "\t0/1");
					break;
				case 2:
					fprintf(f_out, "\t0/0");
					break;
				case NAN_VAL:
					fprintf(f_out, "\t./.");
					break;
			}
		}
		free(record);
		fprintf(f_out, "\n");
	}
	fclose(f_out);
}


void print_help() {
	printf("Help goes here!\n");
}

void print_multiple_flag_error(char* flag) {
	fprintf(stderr, "Error: Multiple instances of %s flag!\n", flag);
}

int main(int argc, char* argv[]) {
	char* in_vcf = NULL;
	char* in_geno = NULL;
	char* in_ind = NULL;
	char* in_snp = NULL;
	bool is_prefix = false;
	
	char* out_vcf = NULL;
	char* out_geno = NULL;
	char* out_ind = NULL;
	char* out_snp = NULL;

	char** sex = NULL; size_t n_sex;
	char** pop = NULL; size_t n_pop;

	static struct option long_options[] = {
		{"help",        no_argument,       NULL, 'h'},
		{"vcf",         required_argument, NULL, 'v'},
		{"geno",        required_argument, NULL, 'g'},
		{"ind",         required_argument, NULL, 'i'},
		{"snp",         required_argument, NULL, 's'},
		{"prefix",      required_argument, NULL, 'p'},
		{"sex",         required_argument, NULL, 'S'},
		{"pop",         required_argument, NULL, 'P'},
		{"out-prefix",  required_argument, NULL, 'o'},
		{"out-vcf",     required_argument, NULL, 'O'}
	};

	size_t alloc_size = 0;  // to store for dynamic string creation;
	while(1) {
		int c = getopt_long(argc, argv, "hp:P:s:i:o:", long_options, NULL);
		if(c == -1) { break; }
		switch(c) {
			case 'h':
				print_help();
				exit(EXIT_SUCCESS);
			case 'v':
				if(in_vcf) {
					print_multiple_flag_error("'--vcf'");
					exit(EXIT_FAILURE);
				}
				in_vcf = optarg;
				break;
			case 'g':
				if(in_geno) {
					print_multiple_flag_error("'--geno' or '--prefix'");
					exit(EXIT_FAILURE);
				}
				in_geno = optarg;
				break;
			case 'i':
				if(in_ind) {
					print_multiple_flag_error("'--ind' or '--prefix'");
					exit(EXIT_FAILURE);
				}
				in_ind = optarg;
				break;
			case 's':
				if(in_snp) {
					print_multiple_flag_error("'--snp' or '--prefix'");
					exit(EXIT_FAILURE);
				}
				in_snp = optarg;
				break;
			case 'p':
				if((in_geno != NULL) || (in_snp != NULL) || (in_ind != NULL)) {
					fprintf(stderr, "ERROR: PACKEDANCESTRYMAP file input information already provided!\n");
					exit(EXIT_FAILURE);
				}
				is_prefix = true;
				alloc_size = strlen(optarg) + PAM_STR_BUF_EXTRA + 1;
				in_geno = (char*)malloc(alloc_size * sizeof(char));
				in_snp = (char*)malloc(alloc_size * sizeof(char));
				in_ind = (char*)malloc(alloc_size * sizeof(char));
				sprintf(in_geno, "%s.geno", optarg);
				sprintf(in_snp, "%s.snp", optarg);
				sprintf(in_ind, "%s.ind", optarg);
				break;
			case 'S':
				if(sex) {
					print_multiple_flag_error("'--sex'");
					exit(EXIT_FAILURE);
				}
				sex = str_split(optarg, ',', &n_sex);
				break;
			case 'P':
				if(pop) {
					print_multiple_flag_error("'--pop'");
					exit(EXIT_FAILURE);
				}
				pop = str_split(optarg, ',', &n_pop);
				break;
			case 'o':
				if((out_geno != NULL) || (out_snp != NULL) || (out_ind != NULL)) {
					fprintf(stderr, "ERROR: PACKEDANCESTRYMAP genotype file output information already provided!\n");
					exit(EXIT_FAILURE);
				}
				alloc_size = strlen(optarg) + PAM_STR_BUF_EXTRA + 1;
				out_geno = (char*)malloc(alloc_size * sizeof(char));
				out_snp = (char*)malloc(alloc_size * sizeof(char));
				out_ind = (char*)malloc(alloc_size * sizeof(char));
				sprintf(out_geno, "%s.geno", optarg);
				sprintf(out_snp, "%s.snp", optarg);
				sprintf(out_ind, "%s.ind", optarg);
				break;
			case 'O':
				if(out_vcf) {
					print_multiple_flag_error("'--out-vcf'");
					exit(EXIT_FAILURE);
				}
				out_vcf = optarg;
				break;
			case '?':
				break;
			default:
				abort();
		}
	}
	// check input integrity
	bool in_pam = (in_geno != NULL) || (in_snp != NULL) || (in_ind != NULL);
	bool out_pam = (out_geno != NULL) || (out_snp != NULL) || (out_ind != NULL);
	if((in_vcf != NULL) && in_pam) {
		fprintf(stderr, "ERROR: both PACKEDANCESTRYMAP and VCF inputs provided!\n");
		exit(EXIT_FAILURE);
	}
	// check output integrity
	if((out_vcf != NULL) && out_pam) {
		fprintf(stderr, "ERROR: both PACKEDANCESTRYMAP and VCF inputs provided!\n");
		exit(EXIT_FAILURE);
	}
	// check input and output alignment
	if((in_vcf != NULL) && (out_vcf != NULL)) {
		 fprintf(stderr, "ERROR: VCF input provided alongside VCF output. No conversion can be performed.\n");
		 exit(EXIT_FAILURE);
	}
	if(in_pam && out_pam) {
		fprintf(stderr, "ERROR: PACKEDANCESTRYMAP input provided alongside PACKEDANCESTRYMAP output. No conversion can be performed.\n");
		exit(EXIT_FAILURE);
	}

	// final run
	if(in_pam) {
		ind_data ind = read_ind_file(in_ind);
		snp_data snp = read_snp_file(in_snp);
		pam_file_reader pfr = pam_file_reader_init(in_geno, &snp, &ind);
		read_pam_header(&pfr);	
		write_vcf(&pfr, out_vcf, &snp, &ind);
		free_ind_data(&ind);
		free_snp_data(&snp);
		close_pam_file_reader(&pfr);
	} else if(in_vcf != NULL) {
		// chekc that population information was provided
		if(pop == NULL) {
			fprintf(stderr, "ERROR: --pop parameter must be provided to convert VCF to PACKEDANCESTRYMAP.\n");
			exit(EXIT_FAILURE);
		}
		if((sex != NULL) && (n_sex != n_pop)) {
			fprintf(stderr, "ERROR: invalid number of '--pop' or '--sex' labels provided!\n");
			exit(EXIT_FAILURE);
		}
		FILE* f_vcf = fopen(in_vcf, "r");
		if(f_vcf == NULL) {
			fprintf(stderr, "ERROR: could not open VCF file.\n");
			exit(EXIT_FAILURE);
		}
		write_snp_and_ind(f_vcf, out_snp, out_ind, pop, sex, n_pop);
		ind_data ind = read_ind_file(out_ind);
		snp_data snp = read_snp_file(out_snp);
		write_geno(f_vcf, &snp, &ind, out_geno);
		free_ind_data(&ind);
		free_snp_data(&snp);
		fclose(f_vcf);
	}
	free(out_geno); free(out_ind); free(out_snp);
	if(is_prefix) { free(in_geno); free(in_ind); free(in_snp); }
	if(sex) {
		for(size_t i = 0; i < n_sex; i++) {
			free(sex[i]);
		}
		free(sex);
	}
	if(pop) {
		for(size_t i = 0; i < n_pop; i++) {
			free(pop[i]);
		}
		free(pop);
	}
}
