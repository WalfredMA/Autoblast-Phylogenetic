# Autoblast-Phylogenetic
generate Phylogenetic tree by cloud blasting

Basically what this program does:

In a folder called result in R's default folder
It generates a loop to:

1. Call on ncbi's blast sites and extract data each turn.
2.Parameters can be set on the top findsimilar means the program will not stop until it got 2 run has same number of hit.
3. From the result file, it extracts data and generate forms in each run's folder
4. It download full sequence fasta file of each aligned  protein.
5. It generate phylogenetic tree from fasta files
6. It substrate species' names we extracted before to replace accession number in the tree. Since the phylotree does not allow same name...I put a index before same species.
7. I put a try loop to prevent the break out when website 404...the max try I set 10 times try...but still will break out sometimes....but we can set to 100 times to prevent this...however...it may result   in system freezing for a very long time...
8. I tried several times...but not too many...so maybe still exist some bugs..
9. I post a sample result here...it took about 1h30mins for 10 terms run
10. If it break out at some step...you keep results generated before and continuo by manually putting i and runtime values. I use name very widely...so you can check which step it 404 by typing name on console..
