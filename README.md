# Consensus

### Librairies
* https://github.com/ksahlin/NGSpeciesID
* https://github.com/asrivathsan/ONTbarcoder

### How to run the sequence_preprocessing.py script

The Anaconda environment used to run the sequence_preprocessing script can be found in the file GenoRobotics_bioinfo_consensus.yml
<br/>Here are the following steps to follow to run it:
 - open the current folder in your terminal
 - type `conda env create -f GenoRobotics_bioinfo_consensus.yml` in the terminal to create the environment GenoRobotics_bioinfo_consensus (Anaconda must be installed)
 - type `conda activate GenoRobotics_bioinfo_consensus` in the terminal to activate the environment
 - type `cd src` in the terminal
 - type `python sequence_preprocessing.py` in the terminal to try to run the main loop of sequence_preprocessing.py
 - to exit the environment, type `conda deactivate` in the terminal

If everything can be run, it means that the environment is correctly installed.
<br/>Well done, you have successfully preprocessed your first dataset ! Now you can try to play with other datasets with different parameters. 👍

### Notes

Currently, a small dataset called `rbcL_Qiagen_tomato_5000.fastq` (containing 5000 sequences, ~5 MB) can be found in the `src` folder. This is a sampled version from the entire dataset `rbcL_Qiagen_tomato.fastq` (containing ~200k sequences, ~204 MB). This dataset is used in the main loop of the sequence_preprocessing.py script in the run example and can be used to check if the environment is correctly installed.

The current version of the sequence_preprocessing.py is a test/incomplete version and everything is subject to change.
