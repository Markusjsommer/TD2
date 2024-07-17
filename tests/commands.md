# Comparing pep outputs
conda activate td
pip install .

TD2.LongOrfs -t ./tests/data/MANE_test_rna.fna --m-start -O ./results/
python ./tests/compare_peps.py ./tests/data/MANE_test_rna.fna.transdecoder_dir/longest_orfs.pep ./results/transcripts.transmark_dir/longest_orfs.pep --diff --matched --internal

# Installing TransDecoder
conda activate td
conda install -c conda-forge -c bioconda transdecoder

TransDecoder.LongOrfs -h
TransDecoder.Predict -h

# Benchmarking time and memory
* memory -> valgrind --tool=massif [program] 
* memory_output -> ms_print massif.out.<pid>
* time -> /usr/bin/time -v [program]

/usr/bin/time -v TransDecoder.LongOrfs -t data/MANE.GRCh38.v1.3.refseq_rna.fna -O results/TD_benchmark/time/ > results/TD_benchmark/time.log
valgrind --tool=massif TransDecoder.LongOrfs -t data/MANE.GRCh38.v1.3.refseq_rna.fna -O results/TD_benchmark/memory/ > results/TD_benchmark/memory.log