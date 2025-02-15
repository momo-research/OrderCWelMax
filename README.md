# OrderCWelMax

## **Directory Structure**
```
OrderCWelMax/
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
│   ├── 
├── makefile                  # makefile for compiling
└── README.md                 # Project documentation
```
---

## **Dependencies**
- **C++ Compiler**: Supports C++17 or higher.
- **Libraries**:
  - Standard Template Library (STL)
  - `<omp.h>`, `<mutex>` et al.

---

## **Build Instructions**
1. Clone the repository:
2. Enter the Source Code directory
3. Compile the project
We provide a Makefile for easy compilation. Simply run the following command:
```
  make
```
You can modify the compilation options in the Makefile if needed.
4. run the program
```
./build/seqic
```
   
