# ADROITrimer
An Automated Pipeline for HIV-1 Envelope Glycoprotein Sequence Design To Achieve Prefusion-Closed Trimers


## Installation

ADROITrimer has been successfully tested on Mac and Linux systems

### Dependencies

- MUSCLE alignment software (https://www.drive5.com/muscle/). Please download and install.
- The R Project for Statistical Computing (https://www.r-project.org/). Please download and install.
  - R library bio3d
    Install by opening R in command line and typing: install.packages("bio3d")


## Run 
ADROITrimer can be run in the command line

### Dependencies
3 files are required to be placed in work directory (please download them from this repository in the required_files directory)
  1. HIV1_FLT_2016_env_PRO_cleanLastCol_final.fasta
  2. 5fyl_trimer.pdb
  3. 5fyl_trimer.rsa

### 3 input arguments are necessary to run the ADROITrimer
  1.  HIV-1 Env target sequence (sequence you would like to be optimized) in FASTA format
  2.  Clade/Subtype of target sequence (Clade can be for instance determined using RIP from LANL: https://www.hiv.lanl.gov/content/sequence/RIP/RIP.html)
  3.  Path to installed MUCLE software (e.g. "/usr/local/bin/muscle")
  
### Execute in the command line
R --no-save < ADROITrmer1.0-1.R target.fasta A /usr/bin/muscle


## Result
Multiple analysis and output files will be provided. The critical DS-SOSIP-RnS stabilized sequences are provided in FASTA format. In particular,
1. With N- and C-terminal tags (N-terminal leader sequence, C-terminal - D7324-tag): output_xxx_DS-SOSIP-RnS_NtCt-tags.fasta (with xxx being a numeric value - time stamp)
2. Without N- and C-terminal tags: output_xxx_DS-SOSIP-RnS (with xxx being a numeric value - time stamp)
  You would have to add signal peptide and any desired C-terminal tag.
  
  
### Example/Test

#### Clade A Strain BG505 testing sequence (please find BG505.W6M.C2.fasta example sequence in test sub-directory) 
R --no-save < ADROITrmer1.0-1.R BG505.W6M.C2.fasta A /usr/bin/muscle


# Contact
Reda Rawi: reda.rawi@nih.gov
