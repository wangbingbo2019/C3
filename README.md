# C3

C3.py runs the C3 method as described in
 
C3: Connect separate Connected Components to form a succinct disease module
by Bingbo Wang, Jie Hu, Chenxing Zhang, Yuanjun Zhou, Liang Yu, Xingli Guo, Lin Gao

Instruction to use the source code:
1. Download the code.
2. Make sure to make the code executable by chmod +x C3.py
3. Run the code (./C3.py in command line or “run C3.py" in ipython).
4. Run the code according the input files description(See below).

# -------------------

Directory Example

contains two input files:
1. seeds_file.txt (list of disease genes) 

2. network_file.txt (human interactome. note that gene IDs should be consistent in the two input files) This network is provided by 'Menche J, Sharma A, Kitsak M, et al. Disease networks. Uncovering disease-disease relationships through the incomplete interactome. Science 2015; 347.'

The following command will generate the C3 genes and connected disease genes and save them in two files)
./C3.py network_file.txt seeds_file.txt
