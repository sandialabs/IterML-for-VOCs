Companion repository to "Combining Machine Learning and Iterative Experiments to Keep Pace with Emerging Viral Variants of Concern"
Thomas Sheffield, Ryan C. Bruneau, Stephen Won, Kenneth L. Sale, Brooke Harmon, Le Thanh Mai Pham. PLOS Computational Biology (2026)



R Packages used:  
Matrix 1.7-2  
data.table 1.17.0  
protr 1.7-5  
parallel 4.4.3  
openxlsx 4.2.8  
Rfast 2.1.5.2  
keras 2.16.0  
numDeriv 2016.8-1.1  
xgboost 1.7.10.1  
ranger 0.17.0  
seqinr 4.2-36  
gplots 3.2.0  
RColorBrewer 1.1-3



Python Libraries:  
python 3.9.15  
keras 2.10.0  
tensorflow 2.10.0  
pandas 2.3.3  
binarymap 0.8  
dms\_variants 1.6.0  



I recommend setting up a virtual environment with the required libraries and selecting this environment as your python interpreter in RStudio global options.
For all functions, source the file it's contained in before running it. I have tried to include key intermediate files as much as possible, including the final
public data NN models used to inform the RF models and select variants. I have also tried to identify where to find important outputs, but not all outputs.



***Data***  


Where to find public data:  
PACE: https://github.com/jbloomlab/SARS-CoV-2-RBD\_DMS/tree/master/results/binding\_Kds/binding\_Kds.csv  
PEx: https://github.com/jbloomlab/SARS-CoV-2-RBD\_DMS/tree/master/results/expression\_meanFs/expression\_meanFs.csv  
PAnti: https://github.com/jbloomlab/SARS-CoV-2-RBD\_MAP\_Crowe\_antibodies/tree/master/results/escape\_scores/scores.csv  



Where to put public data:  
PACE:  Datasets/SARS-CoV-2-RBD\_DMS/results/binding\_Kds/binding\_Kds.csv  
PEx:   Datasets/SARS-CoV-2-RBD\_DMS/results/expression\_meanFs/expression\_meanFs.csv  
PAnti: Datasets/SARS-CoV-2-RBD\_MAP\_Crowe\_antibodies/results/escape\_scores/scores.csv  



Our Data:  
I1: Datasets/Library of 59 antibodies\_ Binding and sequences\_MP 12202022.xlsx  
I2: Datasets/2\_08\_23\_omicron\_ba5\_elisa.xlsx  
Datasets/2\_09\_23\_omicron\_ba5\_elisa.xlsx  
I3: Datasets/Variant selection and binding assay \_MP 09182023.xlsx  
Datasets/Variant selection and binding assay \_MP 08162023.xlsx  



Additional Data:  
Datasets/variant\_sequences.csv  
Datasets/var\_to\_align.fasta (from mai\_seq\_align())  
Datasets/var\_aligned.fasta (from muscle alignment)  
greaney antibodies short.xlsx (from literature)  
output/var\_muts.csv (from var\_muts())  



**Data preprocessing** 
R/ml\_preprocess.R  
R/mai\_antib\_bind.R  
R/conc\_resp.R  
R/combined\_preprocess.R



source("R/ml\_preprocess.R")  
PEx:  ml\_preprocess("express")  
PACE: ml\_preprocess("bind")  
PAnti: ml\_preprocess("log\_escape")



source("R/mai\_antib\_bind.R")  
I1: mai\_antib\_bind("12\_20\_2022")  
I2: mai\_antib\_bind("02\_09\_2023")  
I3: mai\_antib\_bind("09\_18\_2023")  



source("R/combined\_preprocess.R")  
Com\_NN: antibody\_preprocess()  
Com\_Epi: binding\_preprocess()  



Key output locations: input/  
output/mai\_antib\_kds\_\*.csv  



***Public NNs*** 
R/starr\_ml.R  
R/ml\_features.R  
R/ml\_plot.R  
R/ml\_modeler.R  
R/ml\_tuner.R  
R/score\_funs.R  



source("R/starr\_ml.R")  
PEx\_NN:  
Default CV            starr\_ml(datatype = "express", feats = "full", ml = "keras")  
Tuned (Fig. S1C)      starr\_ml(datatype = "express", feats = "full", ml = "keras", ntune = 60)  
Full tuned (Fig. S1D) starr\_ml(datatype = "express", feats = "full", ml = "keras", ntune = 60, tunefold = 5)  
Final Model           starr\_ml(datatype = "express", feats = "full", ml = "keras", nfold = 1, ml.params = list(layer.sizes = c(16, 128), relu.alpha = c(.01, .01)), do.plot = F)



PACE\_NN:
Default CV            starr\_ml(datatype = "bind", feats = "full", ml = "keras")  
Tuned (Fig. S1A)      starr\_ml(datatype = "bind", feats = "full", ml = "keras", ntune = 60)  
Full tuned (Fig. S1B) starr\_ml(datatype = "bind", feats = "full", ml = "keras", ntune = 60, tunefold = 5)  
Final Model           starr\_ml(datatype = "bind", feats = "full", ml = "keras", nfold = 1, ml.params = list(layer.sizes = c(64, 4, 16), relu.alpha = c(.01, .01, .01)), do.plot = F)  



PAnti\_NN:  
Default CV        starr\_ml(datatype = "log\_escape", feats = "full", ml = "keras", comb\_antib = T)  
Tuned (Fig. S1E)  starr\_ml(datatype = "log\_escape", feats = "full", ml = "keras", comb\_antib = T, ntune = 60)  
Full tuned        starr\_ml(datatype = "log\_escape", feats = "full", ml = "keras", comb\_antib = T, ntune = 60, tunefold = 5)  
Final Model       starr\_ml(datatype = "log\_escape", feats = "full", ml = "keras", comb\_antib = T, nfold = 1, ml.params = list(layer.sizes = c(64, 4, 16), relu.alpha = c(.01, .01, .01)), do.plot = F)  



Key output locations: output/ml\_express\*/  
output/ml\_bind\*/  
output/ml\_log\_escape\*/  



***Predictions and variant selection***  
R/ml\_predict.R  
R/site\_chooser.R  



Make multiple predictions:  
source("R/ml\_predict.R")  
Omicron within 2 mutations:   pred\_multi\_driver(inpaths = c("ml\_bind\_keras\_full\_layer.sizes\_c(64, 4, 16)\_relu.alphas\_c(0.01, 0.01, 0.01)\_onefold0/Binding.h5", "ml\_express\_keras\_full\_layer.sizes\_c(16, 128)\_relu.alphas\_c(0.01, 0.01)\_onefold/Expression.h5", "ml\_log\_escape\_keras\_full\_layer.sizes\_c(256, 32, 16)\_relu.alphas\_c(0, 0, 0)\_comb\_onefold/Antibody.h5"), target = "Omicron", mutmax = 2)  
BA4\_BA5 within 2 mutations:   pred\_multi\_driver(inpaths = c("ml\_bind\_keras\_full\_layer.sizes\_c(64, 4, 16)\_relu.alphas\_c(0.01, 0.01, 0.01)\_onefold0/Binding.h5", "ml\_express\_keras\_full\_layer.sizes\_c(16, 128)\_relu.alphas\_c(0.01, 0.01)\_onefold/Expression.h5", "ml\_log\_escape\_keras\_full\_layer.sizes\_c(256, 32, 16)\_relu.alphas\_c(0, 0, 0)\_comb\_onefold/Antibody.h5"), target = "BA4\_BA5", mutmax = 2)  



Select 213 variants for further experiments:  
source("R/site\_chooser.R")  
variant\_selection()  



Key output locations: output/multipreds/  
output/Variant\_Selection\_\*.xlsx  



***Ix\_RF Models***  
R/mai\_antib\_model.R  



source("R/mai\_antib\_model.R")  
I1\_RF CV (Fig 3A) mai\_antib\_model(delta.kd = "wuhan", version = 1, use.col = "blue")  
I2\_RF CV (Fig 3B) mai\_antib\_model(delta.kd = "wuhan", version = 2, use.col = "darkgreen")  
I3\_RF CV (Fig 3C) mai\_antib\_model(delta.kd = "wuhan", version = 3, use.col = "purple4", keep.var.feats = 10)  



Key output locations: output/mai\_antib\_model/



***Com\_NN***  
R/combined\_ml.R  
R/combined\_features.R  
R/combined\_plot.R  
R/ml\_modeler.R  
R/ml\_tuner.R  
R/score\_funs.R  



source("R/combined\_ml.R")  
CV (Fig 5) combined\_ml(ml = "keras")  



Key output locations: output/antibody\*/  



**Com\_Epi**  
R/combined\_ml.R  



source("R/combined\_ml.R")  
CV (Fig. S3)               combined\_ml(datatype = "binding", ml = "epistasis")  
Single Mutation Effects    combined\_ml(datatype = "binding", ml = "epistasis", nfold = 1)



Key output locations: output/binding\_epistasis\*/



***Figs 5 and 6***  
R/site\_chooser.R  
R/binding\_interactions.R



source("R/site\_chooser.R")  
Scaled Sum Plot (Fig 6)    new\_starr\_regions()  

source("R/binding\_interactions.R")  
Region Interaction (Fig 7) binding\_interactions()  



Key output locations: output/binding\_combined\_site\_regions.pdf  
output/binding\_epistasis\_epistasis\_none\_onefold\_deltawuhan\_fixscale/reg\_int\_plot\_final.pdf

