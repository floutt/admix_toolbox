# admix_toolbox
`admix_toolbox` is a collection of programs designed for working with PACKEDANCESTRYMAP and TGENO files. Many of these programs replace (and outperform) functionality traditionally done using `AdmixTools`.

## Installation
These programs can be installed as follows:
```sh
make
```
Users seeking to install this program globally can then simply run:
```sh
make install
```
Users who seek to install the binaries in a local directory can simply change the `PREFIX` parameter. For example, someone who wishes to install these programs in `$HOME/bin/` can run:
```sh
make install PREFIX=$HOME
```

## Programs
### pam_merge
The `mergeit` command is typically used to merge PACKEDANCESTRYMAP files, however, this program loads the contents of entire files into memory. `pam_merge` fixes this issue by reading data on a record-by-record basis.

Information on running the command can be obtained by running `pam_merge -h`:
```
Usage: 'pam_merge [OPTIONS] --prefixes [PREFIX] --out [OUTFILE]' or 'pam_merge [OPTIONS] --pams [GENOFILES] --snps [SNPFILES] --inds [INDFILES] --out [OUTFILE]'
Options:
	-h, --help        Display help message and exit
	-p, --prefixes <list>  Comma-separated list of prefixes of PACKEDANCESTRYMAP files to merge
	-P, --pams <list>      Comma-separated list of PACKEDANCESTRYMAP geno files to merge
	-s, --snps <list>      Comma-separated list of PACKEDANCESTRYMAP SNP files to merge
	-i, --inds <list>      Comma-separated list of PACKEDANCESTRYMAP individual files to merge
	-o, --out <prefix>     Prefix of output PACKEDANCESTRYMAP files
	--ignore-hash          Ignore PACKEDANCESTRYMAP hash check
	--intersect-snps       Return files with the intersection of SNPs from all files, if this flag is not set it will use first file as SNP reference
```

### vcf_pam_conv
`vcf_pam_conv` provides an easy way to convert between VCF and PACKEDANCESTRYMAP files. To convert a PACKEDANCESTRYMAP file set (`example.*`) to a VCF can simply run:
```sh
vcf_pam_conv --geno example.geno --ind example.ind --snp example.snp --out-vcf out.vcf
```
Similarly, to convert a VCF file (`in.vcf`) to a PACKEDANCESTRYMAP file set one can run the following command:
```sh
vcf_pam_conv --vcf in.vcf --sex sex.samples --pop pop.samples --out-prefix out_pam
```
More information on this command can be obtained by running `vcf_pam_conv -h`:
```
Usage: vcf_pam_conv [OPTIONS]
Options:
        -h, --help                      Display help message and exit
        -v, --vcf <filename>            VCF input file
        -g, --geno <filename>           PACKEDANCESTRYMAP genotype file
        -s, --snp <filename>            PACKEDANCESTRYMAP SNP file
        -i, --ind <filename>            PACKEDANCESTRYMAP ind file
        -s, --snp <prefix>              PACKEDANCESTRYMAP file prefix
        -S, --sex <list>                Comma-separated list of sex IDs for VCF file individuals
        -P, --pop <list>                Comma-separated list of population IDs for VCF file individuals
        -o, --out-prefix <prefix>       Prefix for PACKEDANCESTRYMAP output
```

### transpose_geno
`transpose_geno` is a program for converting a TGENO file to the more widely used PACKEDANCESTRYMAP format. While `convertf` (a part of the `AdmixTools` software) can do this task, it requires loading the entire TGENO file into memory. `transpose_geno` allows the user to adjust the amount of memory they seek to use in the conversion process through the `--n-snp` parameter. Higher values use more memory but run faster, lower values use less memory but run slower.

The program can be run as follows:
```sh
transpose_geno --prefix [TGENO PREFIX] --out-geno [OUTPUT PACKEDANCESTRYMAP FILE] --n-snp [NUMBER OF SNP RECORDS TO KEEP IN MEMORY]
```
