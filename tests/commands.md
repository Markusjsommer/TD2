# Comparing pep outputs
conda activate td
pip install .

TD2.LongOrfs -t ./tests/data/MANE_test_rna.fna --m-start -O ./results/
python ./tests/compare_peps.py ./tests/data/MANE_test_rna.fna.transdecoder_dir/longest_orfs.pep ./results/transcripts.transmark_dir/longest_orfs.pep --diff --matched --internal

## Fixed
TD2.LongOrfs -t ./tests/data/MANE_test_rna.fna -O ./results/

# Installing TransDecoder
conda activate td
conda install -c conda-forge -c bioconda transdecoder

TransDecoder.LongOrfs -h
TransDecoder.Predict -h

# Benchmarking time and memory
* memory -> valgrind --tool=massif [program] 
* memory_output -> ms_print massif.out.<pid>
* time -> /usr/bin/time -v [program]

** valgrind massif doesn't seem to generate a massif file ??? need to rely on time data

## TD
conda activate td
/usr/bin/time -v TransDecoder.LongOrfs -t data/MANE.GRCh38.v1.3.refseq_rna.fna -O results/TD_benchmark/time/ > results/TD_benchmark/time.log 2>&1
valgrind --tool=massif --massif-out-file=results/TD_benchmark/massif.out TransDecoder.LongOrfs -t data/MANE.GRCh38.v1.3.refseq_rna.fna -O results/TD_benchmark/memory/ > results/TD_benchmark/memory.log 2>&1
ms_print massif.out

## TD2
conda activate td2
/usr/bin/time -v TD2.LongOrfs -t data/MANE.GRCh38.v1.3.refseq_rna.fna -O results/TD2_benchmark/time/ > results/TD2_benchmark/time.log 2>&1
valgrind --tool=massif --massif-out-file=results/TD2_benchmark/massif.out TD2.LongOrfs -t data/MANE.GRCh38.v1.3.refseq_rna.fna -O results/TD2_benchmark/memory/ > results/TD2_benchmark/memory.log 2>&1
ms_print massif.out

## compare outputs
python ./tests/compare_peps.py ./results/TD_benchmark/time/longest_orfs.pep ./results/TD2_benchmark/time/transcripts.transmark_dir/longest_orfs.pep --diff --matched --internal

## multithreading
TD2.LongOrfs -t data/MANE.GRCh38.v1.3.refseq_rna.fna -@ 32 -O results/multi/