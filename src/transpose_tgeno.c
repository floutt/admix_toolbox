#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <getopt.h>
#include <dotgeno.h>

#define PAM_STR_BUF_EXTRA 6
#define MIN(a, b) ((a) < (b) ? (a) : (b))

// loading stuff
void print_progress(char* message, size_t cur, size_t max, size_t len) {
	float prop = (1.0 * cur)/max;
	float percent = prop * 100;
	size_t num_fill = (size_t)roundf(prop * len);
	size_t num_empty = len - num_fill;
	printf("%s: [", message);
	for(size_t i = 0; i <= num_fill; i++) {
		printf("=");
	}
	for(size_t i = 0; i < num_empty; i++) {
		if(i == 0) { printf(">"); continue; }
		printf(".");
	}
	printf("]%.0f%\r", percent);
	fflush(stdout);
	if(cur == max) { printf("\n"); }
}

int main(int argc, char* argv[]) {
	bool intersect_snps = false;
	char* ind_file = NULL;
	char* snp_file = NULL;
	char* geno_file = NULL;
	char* out_geno = NULL;
	size_t n_snp = 0;
	static struct option long_options[] = {
		{"prefix",       required_argument, NULL, 'p'},
		{"out-geno",       required_argument, NULL, 'o'},
		{"n-snp",          required_argument, NULL, 'n'},
		{0,                0,                 0,      0}
	};

	while(1) {
		int c = getopt_long(argc, argv, "hp:P:s:i:o:", long_options, NULL);
		if(c == -1) { break; }
		switch(c) {
			case 'p':
				if(geno_file) {
					fprintf(stderr, "ERROR: PACKEDANCESTRYMAP files already provided!\n");
					exit(EXIT_FAILURE);
				}
				snp_file = (char*)malloc(sizeof(char) * (strlen(optarg) + PAM_STR_BUF_EXTRA + 1));
				ind_file = (char*)malloc(sizeof(char) * (strlen(optarg) + PAM_STR_BUF_EXTRA + 1));
				geno_file = (char*)malloc(sizeof(char) * (strlen(optarg) + PAM_STR_BUF_EXTRA + 1));
				sprintf(snp_file, "%s.snp", optarg);
				sprintf(ind_file, "%s.ind", optarg);
				sprintf(geno_file, "%s.geno", optarg);
				break;
			case 'o':
				if(out_geno) {
					fprintf(stderr, "ERROR: Output file already provided!\n");
					exit(EXIT_FAILURE);
				}
				out_geno = optarg;
				break;
			case 'n':
				if(n_snp) {
					fprintf(stderr, "--n-snp value already provided!\n");
					exit(EXIT_FAILURE);
				}
				n_snp = (size_t)atoi(optarg);
			case '?':
				break;
			default:
				abort();
		}
	}

	if(n_snp == 0) { n_snp = 1; }

	snp_data snp = read_snp_file(snp_file);
	ind_data ind = read_ind_file(ind_file);
	tgn_file_reader tfr = tgn_file_reader_init(geno_file, &snp, &ind);
	read_tgn_header(&tfr);
	pam_file_writer pfw = pam_file_writer_init(out_geno, &snp, &ind);
	write_pam_header(&pfw, &snp, &ind);
	uint8_t* record;
	size_t max_i = (size_t)ceil((1.0 * snp.length) / n_snp);
	print_progress("Converting", 0, max_i, 20);
	for(size_t i = 0; i < max_i; i++) {
		print_progress("Converting", i+1, max_i, 20);
		size_t offset = i * n_snp;
		size_t n = MIN(n_snp, snp.length - offset);
		uint8_t** snp_mat = (uint8_t**)malloc(n * sizeof(uint8_t*));
		for(size_t j = 0; j < n; j++) { snp_mat[j] = (uint8_t*)malloc(ind.length * sizeof(uint8_t)); }
		size_t ind_i = 0;
		while(record = read_tgn_record(&tfr)) {
			for(size_t j = 0; j < n; j++) {
				snp_mat[j][ind_i] = record[j + offset];
			}
			ind_i++;
			free(record);
		}
		for(size_t j = 0; j < n; j++) { write_pam_record(&pfw, snp_mat[j]); free(snp_mat[j]); }
		free(snp_mat);
		goto_ind_tgn(&tfr, &ind, ind.ind_id[0], ind.population[0]);
	}
	free_snp_data(&snp);
	free_ind_data(&ind);
	close_tgn_file_reader(&tfr);
	close_pam_file_writer(&pfw);
	free(snp_file); free(ind_file); free(geno_file);
	printf("Done!\n");
}
