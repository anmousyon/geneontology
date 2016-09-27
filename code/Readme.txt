1. Install following packages:

	a.Install following through Bioconductor
	 - Here is a link that explains how to install packages through Bioconductor
		GO.db
		GOSim
		GOSemSim
		GOstats
		RBGL

	b. Install normally. In RStudio, go to Tools -> Install Packages
		arules
		
2. From within RStudio, set the working directory to the "code" folder.
	For example, in Mac OS X if location of code folder is "~/GO_Asssociation_Paper/code", then do the following:
	setwd("~/GO_Asssociation_Paper/code") 
	In Windows, you might have to use use double slashes:
	setwd("C://temp//....") <- I am not 100% sure, as I don't use Windows

3. Go to code folder and run the findAssociation.R file by sourcing it.
	source("findAssociation.R")
	
4. Be patient, it can take upto 30 minutes. If you want a quick check, reduce size of the gene_association/gene_association.sgd.csv file.

5. YOUR TASK:
	We need to replace gene_association.sgd.csv file with human annotation file. 
	The list of annotations is at: http://geneontology.org/page/download-annotations
	Direct link for human annotations: http://geneontology.org/gene-associations/goa_human.gaf.gz