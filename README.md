# 🔄OrderCWelMax

## **📑Directory Structure**
```
Source Code/
├── Include/                  # Header files
│   ├── anyoption.h           # Command-line option parsing
│   ├── graph.h               # Graph data structure
│   ├── head.h                # Common headers and macros
│   ├── ingraph.h             # Reverse Influence Sampling
│   ├── mc.h                  # Monte Carlo simulation
│   ├── memoryusage.h         # Memory usage tracking
│   ├── random_utils.h        # Random number generator
│   ├── timgraph.h            # RIS algorithms
├── Src/                      # Source files
│   ├── anyoption.cpp         # Implementation for command-line parsing
│   ├── comic.cpp             # Main program
│   ├── ingraph.cpp           # Implementation for RIS algorithms
│   ├── mc.cpp                # Implementation for Monte Carlo simulation
│   ├── random_utils.cpp      # Implementation for random number generation
│   ├── timgraph.cpp          # Implementation for RIS 
├── Dataset/                  # Datasets         
│   ├── lastfm_wc.txt         # LastFM dataset   n=2100   m=12717
│   ├── flixster_t01.txt      # Flixster dataset n=29358  m=212614
│   ├── attribute.txt         # Dataset features
├── config.txt                # Project configuration
├── makefile                  # Makefile for compiling
└── README.md                 # Project documentation
```
---

## **🧲Dependencies**
- **C++ Compiler**: Supports C++17 or higher.
- **Libraries**:
  - Standard Template Library (STL)
  - `<omp.h>`, `<mutex>` et al.

---

## **🚀Build Instructions**
1. Enter the "Source Code" directory
2. Compile the project. We provide a Makefile for easy compilation. Simply run the following command:
```
  make
```
3. You can modify the compilation options in the Makefile if needed.
4. Run the program
```
./build/seqic
```
## **⚙Configuration**
Project configuration is in config.txt
```
algo : 5             # Choose algorithms
phase : 1            # Phase in some algorithms
ignoreBSeeds : 0     # Whether to read B seeds or not
selectBSeeds : 0     # Whether to randomly select B seeds or not
epsilon : 0.5        # Trade-off between efficiency and solution quality (smaller values improve quality but increase runtime)
theta : 100000       # Default number of RR sets/trees to generate
kA : 10              # Size of A seeds set
kB : 10              # Size of B seeds set
runs : 1             # Number of runs
qao : 0.90           # Global Adoption Parameters
qbo : 0.90
qab : 0.20
qba : 0.60
```

### **Algorithms**
Set option "algo" in config.txt :
```
1 - compute total adoption count with sequential propagation.
2 - compute node count with sequential propagation.
3 - mine seeds using IMM-q algorithm
4 - mine seeds using IMM-NoCompe
5 - estimate SA factor (node count)
6 - estimate SA factor (adoption count)
7 - mine seeds using OPPRT (A frist)
8 - mine seeds using OPPRT (B frist)
9 - mine seeds using Greedy (A first)
10 - mine seeds using Greedy (B first)
11 - mine seeds using High degree
12 - mine seeds using PageRank
```
### **Datasets**
Please modify "attribute.txt" accordingly after your change the dataset.
