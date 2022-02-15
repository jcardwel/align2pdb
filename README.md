# align2pdb
This script will take a multiple protein alignment and annotate it with secondary structure data.  It uses ENSEMBL IDs to match the sequence to the pdb IDs. The protein alignment sequences should be labeled with ensembl IDs.  To do this, it queries and parses results from PBDe, RCSB PDB, and SIFTS.  

# AWS vs. Standalone
There are two versions of the script:
1. AWS - this script is the current version running on the AWS-powered website.  It is designed to run on an uploaded fasta through the Flask website and output to the website.  Downloaded files are all behind the scenes on the AWS instance.
2. Standalone - this is downloaded and used on an independent machine.  Given an input directory and chain file, it will scan subdirectories for what it will assume to be fasta-formatted input files.  It will create, download, and output files to the current working directory, and write permissions will be needed.  Compiled with PyInstaller.  It is executed by invocing the executable directly or by adding the align2pdb/Standalone/alignment_to_dssp_standalone directory to your $PATH.  Example usage: _alignment_to_dssp_standalone -d <fasta_input_directory> -c <chain_file>_

# Website
URL to the AWS-powered Flask website:
http://54.147.54.166:8081


