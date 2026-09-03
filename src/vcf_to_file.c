#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

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
	// COMPLETE
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
