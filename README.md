Command Line Arguments 
- Arg 1 (tp53.txt) : text file including protein query sequence of any length in FASTA format. 
- Arg 2 (gt-tp53.txt) : text file including secondary structure annotation of corresponding protein.

Requirements 
- numpy==1.21.4

1. If you have no ground truth data and/or if you just want to predict the sequence of the second structural elements of a given amino acid sequence, then you can use this option:<br>
`python chou-fasman.py tp53.txt`

Output: <br>
An output.txt file will be created, which includes header of the protein; helix, strand, turn assignments (initial extended hits) and final prediction.


2. If you have the ground truth data of the sequence of the second structural elements of a given amino acid sequence, you can provide this data as a second command line argument to view the test results:<br>
`python chou-fasman.py tp53.txt gt-tp53.txt`

Output:<br>
An output.txt file will be created, which includes header of the protein; helix, strand, turn assignments (initial extended hits) and final prediction. Also confusion matrix; precision, recall, accuracy, and F1-score metrics, for each SS element and overall accuracy score will be displayed on the print screen.


<p align="center">
  <img src="/images/out.PNG">
</p>
