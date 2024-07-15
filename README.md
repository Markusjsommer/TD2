# TD2
A tool to find protein coding ORFs

NOTE: this readme is heavily based on the TransDecoder wiki: https://github.com/TransDecoder/TransDecoder/wiki


TD2 identifies likely coding sequences based on the following criteria:

- a minimum length open reading frame (ORF) is found in a transcript sequence
- a likelihood score computed by PSAURON is > 0.5
- the above coding score is greatest when the ORF is scored in the 1st reading frame as compared to scores in the other 2 forward reading frames and the 3 reverse frames.
TODO ask brian if there is a reason to not use reverse frames here as well. should reduce false positives of short orfs within real orfs.
- if a candidate ORF is found fully encapsulated by the coordinates of another candidate ORF, the longer one is reported. However, a single transcript can report multiple ORFs (allowing for operons, chimeras, etc).




Including homology searches as ORF retention criteria
We suggest using MMseqs2 for increased speed and sensitivity in protein search.
An example command would be like so: 

mmseqs databases UniProtKB/Swiss-Prot swissprot tmp
mmseqs easy-search transdecoder_dir/longest_orfs.pep swissprot alnRes.m8 tmp -s 7.0

The sensitivity of this search can be adjusted with -s. We suggest a setting of -s 7.0 for a very sensitive search.

MMseqs2 can be used to search against a wide variety of databases:
Usage: mmseqs databases <name> <o:sequenceDB> <tmpDir> [options]

  Name                	Type      	Taxonomy	Url
- UniRef100           	Aminoacid 	     yes	https://www.uniprot.org/help/uniref
- UniRef90            	Aminoacid 	     yes	https://www.uniprot.org/help/uniref
- UniRef50            	Aminoacid 	     yes	https://www.uniprot.org/help/uniref
- UniProtKB           	Aminoacid 	     yes	https://www.uniprot.org/help/uniprotkb
- UniProtKB/TrEMBL    	Aminoacid 	     yes	https://www.uniprot.org/help/uniprotkb
- UniProtKB/Swiss-Prot	Aminoacid 	     yes	https://uniprot.org
- NR                  	Aminoacid 	     yes	https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
- NT                  	Nucleotide	       -	https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA
- GTDB                  Aminoacid	     yes	https://gtdb.ecogenomic.org
- PDB                 	Aminoacid 	       -	https://www.rcsb.org
- PDB70               	Profile   	       -	https://github.com/soedinglab/hh-suite
- Pfam-A.full         	Profile   	       -	https://pfam.xfam.org
- Pfam-A.seed         	Profile   	       -	https://pfam.xfam.org
- Pfam-B              	Profile   	       -	https://xfam.wordpress.com/2020/06/30/a-new-pfam-b-is-released
- CDD                   Profile                -        https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml
- eggNOG              	Profile   	       -	http://eggnog5.embl.de
- VOGDB                 Profile                -        https://vogdb.org
- dbCAN2              	Profile   	       -	http://bcb.unl.edu/dbCAN2
- SILVA                 Nucleotide           yes        https://www.arb-silva.de
- Resfinder           	Nucleotide	       -	https://cge.cbs.dtu.dk/services/ResFinder
- Kalamari            	Nucleotide	     yes	https://github.com/lskatz/Kalamari

Additionally, MMseqs2 can be run with on multiple servers using MPI. Detailed notes for MMseqs2 usage can be found here: https://github.com/soedinglab/mmseqs2/wiki