# ExactL3_PPI
The ExactL3 formulations are some better normalized link predictions for Protein-Protein Interaction (PPI) Networks.

This GitHub repo is for an academic conference submission, it includes the introduced algorithms, the experiments to generate all the results data, a jupyter notebook that formats all the result images, and a command-line script to run ExactL3_1, ExactL3_2, and the original L3 link predictor (Kovács, 2019).

# Requirement
Language: ```Python```. Python Libraries (to run the program): ```numpy```, ```pandas```

# Usage
Workflow: example_PPI.txt (input file) => ExactL3_cmd.py => example_PPI_ExactL3_1.txt (output file)

To run ExactL3_cmd.py in the terminal (assume Windows, cmd):
```python ExactL3_cmd.py {input file path} {output file path} {link predictor} {number of CPU core}```

Example to run ExactL3_1:
```python ExactL3_cmd.py ./example_PPI.txt ./example_PPI_ExactL3_1.txt ExactL3_1 1```

Supported Link Predictors: ```ExactL3_1```, ```ExactL3_2```, ```L3```, ```CN```

For examples to work with our Python script, see ```./example.py```. Documentations are included as comments in the script.

# Misc
The data in our paper is generated using the script ```./src/notebook/dataGen.py```, and the images in our paper are generated based on our generated data using the jupyter notebook ```./src/notebook/ppiLPred.ipynb```. Note that except the parsed datasets in ```./src/data/parsed/```, no data is included since it is too large (roughly > 100 GB).

# Docs
For more explanation how our python functions realize the algorithm, see [here](docs/docs.md)
