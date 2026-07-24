Tools in this folder support the W-SLDA infrastructure.  
To compile all tools, execute
```bash
make
```
The tools will be compiled to the subfolder `./bin`.

# List tools
* `td-memory.py`: script that estimates the amount of memory and the minimal number of GPUs needed to run given TD calculations. 
* `high-frequency-filter.py`: simple script showing filter function implemented in TD codes and their impact on the example signal.
* `wdata2checkpoint.c`: Converts W-DATA results into a checkpoint file that can be used as a starting point for the self-consistent process.
* `dpca2wdata-hard.c`:  converts dpca file into wdata set (only for legacy mode).
* `dpca2wdata.c`: creates wtxt metadata file for existing dpca files (only for legacy mode).

