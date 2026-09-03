#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <inttypes.h>
#include <dotgeno.h>

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
	if((n_elems - 9) != ind_len) {
		fprintf(stderr, "ERROR: invalid number of individuals provided. Expected %zu, got %zu\n", n_elems - 9, ind_len);
		exit(EXIT_FAILURE);
	}

	FILE* f_ind = fopen(ind_file, "w+");
	for(size_t i = 9; i < n_elems; i++) {
		char* sex = "U";
		if(ind_sex) { sex = ind_sex[i-9]; }
		fprintf(f_ind, "%s\t%s\t%s\n", elems[i], sex, ind_pop[i-9]);
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
		if((n_elems - 9) != ind->length) {
			fprintf(stderr, "ERROR: invalid number of individuals provided in VCF reecord. Expected %zu, got %zu\n", n_elems - 9, ind->length);
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
			record[i - 9] = dosage;

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
