SOUPE package
================

<!-- badges: start -->

<!-- badges: end -->

The package SOUPE (Sample Optimisation Using Phylogeny and Ecology) aims at
testing and optimizing the phylogenetic and ecological representativeness of a
taxonomical sample in the context of evolutionary biology studies. The method
uses Faith's Phylogenetic Diversity index (PD_Faith) to assess the quality of a
sample. PD_Faith represents the quantity of evolutive history contained in a
sample (Faith 1992). Based on a phylogenetic tree, it calculates the sum of
branch length needed to link all the taxa of the sample. PD_Faith is a relative
value : a higher PD_Faith indicates a more phylogenetically diverse sample. 
Therefore, a sample representative of a studied group, or master-group, has to
maximise PD_Faith compared to other samples of the same size. 

Our method present two possible implementations. The a posteriori process
assesses the quality of an existing sample while the a priori process helps
to create a representative sample according to the phylogeny and ecological
parameters. By "ecological", we mean all non-phylogenetic traits (e.g.,
ecological, phenotypic etc).

In this tutorial, we present both implementation step by step.

## Installation

The first step is to install the package. To do so, download the zip file of the
package on our github page. Do not unzip the file, and directly go to R to
install it in the package installer interface, by selecting the option “package
archive file”, then browsing to the location of the zip file.

Once the package is installed, activate the library :

``` r
library(SOUPE)

```

## A posteriori process : assessing the quality of an existing sample.

PD_Faith is a relative value. In the a posteriori implementation, the tested
sample is compared to a number of random subsets of equal size, generated from
the phylogeny of the master-group. The tested sample is of good quality if its
PD_Faith is higher than the mean of the PD_Faith of the random samples.
If the tested sample was created under limited availability of taxa, it is
possible to give in input the list of available taxa. The random samples will
be created within this list, to ensure that their PD_Faith are comparable.

This process can also evaluate the ecological coverage of the tested
sample, by calculating the number and proportions of species in each traits'
characteristics. It's also possible to assess the ecological representativeness
of the tested subset, i.e. testing if each ecological categories is well 
distributed on the phylogeny. The complete code for this approach is given in
the Ecological traits section.


### Import your files

The files needed to run the code are :

- the phylogeny of the master group, with branch lengths

``` r
plot(ex_tree)
```


- the list of taxa in the tested sample

``` r
print(ex_tested_sample)
```

- if needed, a list of available taxa, that will be used in the creation of
  random samples, so that they can be compared.

``` r
print(ex_available_taxa)
```

- for the ecological coverage, a traits table. It should be composed
  of a column “Species”, and give for each one the discrete traits
  characteristics to be considered.

``` r
print(ex_traits_table)

```

**Advice:**

**- When importing your files, be careful to make sure of the consistency of the
  taxa names between the tested sample, the available taxa, the traits table,
  and the phylogeny. To check if the names match, use `setdiff` function.**

``` r
setdiff(ex_tested_sample, ex_tree$tip.label) 
#returns the elements of the tested_sample that are not in the phylogeny.

```


### Phylogenetic representativeness

Firstly, run the `compost` function. Be careful to store its output in a 
variable. This function calculates PD_Faith of the tested sample (given
in the `sample` argument of the function) and of a number (argument `n`) of
randomly generated samples of the same size, based on the phylogeny of the
master-group (argument `tree`). The random sample can be generated from the
entire master-group (write "tree$tip.label" in the `list` argument), or from a
list of available taxa given by the user (directly given in the `list` argument).
Make sure to store the result of the function in a variable, to use it in the 
next step.


``` r
output_compost <- compost(tree = ex_tree, 
                          sample = ex_tested_sample,
                          list = ex_tree$tip.label, 
                          n = 20)
```

**Advice:**

**- How many sample to generate? Firstly, calculate the number of possible
  random samples that can be generated with the following formula.**
  
**N = (L!)/(T!(L-T)!)**

**In R, it can be calculated with the function `choose(T, L)`.**

**N is the number of samples possible, L is the size of the tested sample 
  (i.e., number of taxa), T is the number of available taxa (i.e., either
  number of taxa in the master-group (tree$tip.label), or number of taxa in the
  list of available taxa).**

**A threshold of 1,000 was set. If N is lower than 1,000, generate the maximal
  number of possible random - 1 (to account for the fact that one of them is the
  sample you are testing). Beyond 1,000 possible samples, it is not possible to
  generate all the random samples for the comparison (too long computation time).
  The number of sample to generate can be 1,000 or more, depending on the size
  of the phylogeny and the tested sample. It should be high enough for the
  comparison to make sense. Here, N = 21, so that we choose to generate 21-1=20
  samples. All the possible samples of this size will therefore be compared.**


Then, use the function `result.compost` to analyse the distribution of PD_Faith
values. In input, give the phylogeny (`tree` argument) and the output of the
compost function, that has been stored in a variable at the previous step (`data`
argument). The result can be visualized with a boxplot (`boxplot` argument), and
a one-sample one-tailed t-test can be conducted to check if the tested sample is 
significantly higher than the mean (`ttes` argument). It is possible to retrieve
the best random sample (`best_random_sample` argument), and to save it in a csv
file if a save location is given in the `save_best_random_sample` argument of
the function.

``` r
result <- result.compost(tree = ex_tree,
                        data = output_compost,
                        boxplot = TRUE,
                        ttest = FALSE,
                        best_random_sample = TRUE)
```

Three points appear on the boxplot : the black point is the mean, the red point 
is the best PD_Faith value, and the blue point is the PD_Faith of the
tested sample. In this example, the tested sample is below the mean. The tested
sample is not of good phylogenetic quality. If the tested sample was above the 
mean, and the t-test significant, then the tested sample would have been of 
good phylogenetic quality. Here, no t-test was conducted, because all possible
samples were generated at the previous step by the function compost.

A non-representative sample can be improved by comparing it to the best random
sample (red point), or to a sample created with the a priori process.


### Ecological coverage

In addition to its phylogenetic quality, it is possible to have information on
the ecological coverage of the tested sample. To do so, use the `metrics`
function. It calculates the proportions and number of species of the tested
sample (`sample` argument) in each traits, from a traits table given in input
(`traits_table` argument). 


``` r
# Traits table
print(ex_traits_table)

# Ecological representativity
test_metrics <- metrics(tree = ex_tree,
                        sample = ex_tested_sample,
                        faith = TRUE,
                        percentage = TRUE,
                        traits = TRUE,
                        traits_table = ex_traits_table, 
                        sub_phylo = TRUE)

```

**Advice:**

**- The function `metrics` can also be used to calculate PD_Faith (`faith`
  argument) and PPR (Percentage of Phylogeny Represented; `percentage` argument)
  of a sample based on the phylogeny of the master-group (`tree` argument), or
  to retrieve and save its phylogeny, pruned from the phylogeny of the
  master-group (`sub_phylo` argument). The phylogeny can be saved as a nexus
  file if a save location is given in the `save_sub_phylo` argument of the
  function.**


#### Ecological representativeness

In addition to the ecological coverage of a sample, one might be interested in
its ecological representativeness. While the `metrics` function gives the 
count and proportion of species in each categories of each ecological traits
present in a traits table, the following code informs about whether the species
that were chosen in each ecological traits are well distributed on the phylogeny, 
and representative of the distribution of the trait in the master-group. It uses
the functions of the a posteriori process, with some data transformation.

For this code, the needed files are : the phylogeny of the master-group, the 
tested sample, and a traits table in the same format as for the `metrics` 
function. Be careful to the consistency of the taxa names between the tested
sample, the traits table, and the phylogeny, and check the absence of NA. 

```
#needed files
test_tree #master-group tree
test_sample #tested sample
test_traits #traits table

#list with all the ecological traits of the traits_table 
all_param <- colnames(test_traits)
all_param

#choose the ecological traits of interest, for example here we chosed "Habitat"
#/!\ only one ecological trait can be chosen at once
n <- which(all_param == "Habitat")

#retrieve the different categories of the chosen ecological traits
categories <- test_traits[[n]] %>% unique()
categories <- unlist(categories)
categories <- categories[! categories %in% NA]

#find the category of each species in the tested sample for this traits
param_tested_sample <- test_traits[[n]][match(test_sample, test_traits$Species)]
names(param_tested_sample) <- test_traits$Species[match(test_sample,
                                                        test_traits$Species)]

#divide the tested sample by ecological category to create one sub-sample
#per category

sub_sample <- list()
for (i in 1:length(categories)) {
  sub_sample[[i]] <- names(param_tested_sample)
                          [param_tested_sample == categories[[i]]]
}
names(sub_sample) <- categories

#divide the master-group species list by ecological category to create one
#sub-sample per category

sub_tree <- list()
for (i in 1:length(categories)) {
  sub_tree[[i]] <- test_traits$Species[test_traits[[n]]==categories[[i]]]
}
names(sub_tree) <- categories

#check the phylogenetical representativeness of the sample for each ecological
#category with the functions of the a posteriori process
length(categories) #gives the number of boxplot that will be created

#optional, split R plot window to display several plots side by side
par(mfrow = c(1,3)) 

#loop
for (i in 1:length(categories)) { #for each category
  
  #retrieve the length of the corresponding subtree
  length_sub_sample_i <- length(sub_sample[[i]]) 
  
  #retrieve the length of the corresponding subsample
  length_sub_tree_i <- length(sub_tree[[i]])
  
  #how many random samples can be generated
  possible_nb <- choose(length_sub_tree_i, length_sub_sample_i) 
  
  if (possible_nb <= 1000) {
    nb_random_samples = possible_nb-1
  }
  if (possible_nb > 1000) {
    nb_random_samples = 1000
  }
  output_compost <- compost(test_tree,
                            sample = sub_sample[[i]],
                            list = sub_tree[[i]],
                            n = nb_random_samples) #compost function
                            
  result.compost(tree = test_tree,
                 data = output_compost,
                 boxplot = TRUE,
                 ttest = TRUE,
                 best_random_sample = FALSE) #result.compost function
                 
  title(sub = categories[[i]])
}

par(mfrow = c(1,1)) #put R plot window back to normal parameter
```


With the a posteriori process, a sample can therefore be easily evaluated for
its phylogenetic and ecological representativeness and coverage. 



## A priori process : What is the best sample for a study?

While the a posteriori process focused on the evaluation of an existing sample, 
the second implementation, namely the a priori process, is an help to the 
creation of a relevant sample for a study.

A phylogenetically relevant sample is a sample that maximises PD_Faith 
for a given number of taxa. The a priori process builds samples iteratively in
order to maximise their Faith indices relatively to the phylogeny of the
master-group. At each step of the iteration, the taxa maximising the increase of
the PD_Faith value of the sample is chosen. As for the a posteriori process, it
is possible to restrict the availability of taxa. Ecological traits can also be
taken into account, resulting in a sample that is optimised for both the
phylogeny and the ecology.


### Import your files

The files needed to run the code are :

- the phylogeny of the master group

``` r
plot(ex_tree)
```

- if needed, the list of available taxa 

``` r
print(ex_available_taxa)

```

- if needed, the traits_table for ecological optimisation. 
``` r
print(ex_traits_table)

```

**Advice:**

**- When importing your files, be careful to make sure of the consistency of the
  taxa names between the tested sample, the available taxa, the traits table,
  and the phylogeny. To check if the names match, use `setdiff` function.**

``` r
setdiff(ex_tested_sample, ex_tree$tip.label) 
#returns the elements of the tested_sample that are not in the phylogeny.

```

### Phylogenetic optimisation

The function `bestsample` is used for phylogenetic optimisation. It allows the
creation of a sample with a maximal PD_Faith for the chosen number of species
(`type = "number"` and `n` arguments), based on the master-group phylogeny
(`tree` argument). As for the a priori process, it’s possible to limit taxa
availability (`list` argument).

**Advices :**

**- To run this code, a taxa has to be chosen to start the iteration (`init` 
  argument of the function). We advise to chose the most “basal” one, or the one
  with the longest branch length.**
  
**- If key species have to be included in the sample, they can be chosen as the
  starting point of the iteration (`init` argument).**
  
**- Sometimes, several taxa equally increase PD_Faith of the sample. The
  algorithm chooses one of them and can store the other ones if a history is 
  requested (`history` argument). It can allow the modification of the created
  sample if needed.**

**- We advise to store the output of the function in a variable to easily access
  the sample composition, its PD_Faith and PPR, and the history, using $. The
  phylogeny of the created sample can also be retrieved (`sub_phylo` argument),
  and saved as a Nexus file when a file path is given in the argument
  `save_sub_phylo`.** 


``` r
output_best_sample <- bestsample(tree = ex_tree,
                                 init = c("G"),
                                 type = "number",
                                 n = 5,
                                 list = ex_available_taxa,
                                 history = TRUE,
                                 sub_phylo = TRUE)
                          
output_best_sample$Best_sample

```

PD_Faith can be expressed as a number but also as a percentage. This
percentage corresponds to the quantity of phylogeny represented in a sample, and
is  calculated as the ratio between PD_Faith of the sample and the total
length of the master-group phylogenetic tree. It is referred hereafter as the
PPR (Percentage of Phylogeny Represented).
With the function bestsample, it is possible to create an optimised sample that
reaches a PPR (`type = "percentage"` and `n` arguments) instead of a number of
taxa.

``` r
output_best_sample2 <- bestsample(tree = ex_tree,
                                  init = c("G"),
                                  type = "percentage",
                                  n = 75,
                                  list = ex_available_taxa,
                                  history = TRUE)
                          
#Here, the function creates a sample that represents 75% of the phylogeny.

output_best_sample2$Best_sample

```


### Ecological optimisation

The function `opti.traits` is used to create a sample optimised for both the 
phylogeny and ecological traits. 
To do so, a traits table is needed. It should be composed of a column “Species”,
and give for each species discrete ecological characteristics to be included.

``` r
print(ex_traits_table)
```

The sample is built in order to maximise PD_Faith while approaching requested
proportions in each ecological trait.

The first step is to chose the proportions to reach for each category of each 
trait.

``` r
proportions <- list(Habitat = c("Sea" = 0.5,
                                "Forest" = 0.5,
                                "Mountain" = 0),
                    Feeding.style = c("Carnivorous" = 0.33,
                                      "Omnivorous" = 0.33,
                                      "Herbivorous" = 0.34),
                    Size = c("Big" = 0.33,
                             "Medium" = 0.33,
                             "Small" = 0.34))
```

**Advices :**

**- When choosing the proportions, make sure that they are realistic compared to
  what is found in the master-group or in the available taxa.**

**- Keep in mind that when several parameters are requested, the created sample
  will most likely only approach the requested proportions, because a compromise
  is done between the ecological parameters and the phylogeny.**


The second step is to weight the traits against each other. It allows to set 
which traits are more important to optimise.

``` r
traits_weight <- c(0.3, 0.3, 0.3)
#Here, each traits has the same weight
```

The third step is to weight the ecology against the phylogeny with the parameter
delta. This parameter quantify whether the sample is optimised preferentially
for the ecology or the phylogeny. Its value is comprised between 0 and 1 :
0 means that only the phylogeny is taken into account, giving the same result
as the `bestsample` function; 1 means that only the ecology is taken into
account, and the sample is not optimised for the phylogeny.

``` r
delta <- 0.5
#Here, the ecology and the phylogeny are equally taken into account.
```

**Advice:**

**- When choosing the delta parameter, keep in mind that if it is set toward
  more ecology, the sample will be less optimised according to the phylogeny. 
  The created sample will be representative of the ecology, but less of the
  phylogeny.**
  
Once the proportions to reach and the relative weights are chosen, the function
`opti.traits` can be run.

``` r
output_optitrait <- opti.traits(tree = ex_tree,
                                n = 5,
                                list = ex_available_taxa,
                                keyspe = c("G"),
                                traits_table = ex_traits_table,
                                proportions = proportions,
                                traits_weight = traits_weight, 
                                delta = delta, 
                                history = FALSE, 
                                sub_phylo = TRUE)

```

**Advice:**

**- As for the bestsample function, an history is available if requested.**

**- Key species can be included to the sample with the `keyspe` option.**

**- Store the output of the function in a variable to easily access the sample
  composition, its PD_Faith, PPR and ecological proportions, and the history,
  using $. The phylogeny of the created sample can also be retrieved 
  (`sub_phylo` argument), and saved as a Nexus file if a file path is given in
  the argument `save_sub_phylo`.**

  
Once the sample has been created with the function `bestsample` or `opti.traits`, 
we advise to check its relevance for the study. If needed, the `history` options
can be used to modify the sample. After any modification, the sample can be
tested again with the a posteriori process (`compost`, `result.compost`, and
`metrics` functions) to check its phylogenetic and ecological representativeness.

``` r
output_compost_bestsample <- compost(tree = ex_tree, 
                                     sample = output_best_sample$Best_sample,
                                     list = ex_available_taxa, 
                                     n = 20)
                                     
result_bestsample <- result.compost(tree = ex_tree, 
                                    data = output_compost_bestsample, 
                                    boxplot = TRUE, 
                                    ttest = TRUE, 
                                    best_random_sample = FALSE)

```
Here, compared to the sample that was tested in the example for the a 
posteriori process (`ex_tested_sample`), the created sample is a good
representative of the phylogeny of the master-group. The initially tested
sample can be compared to the sample created with `bestsample` to improve its
phylogenetic representativeness.


## Visualisation of SOUPE results

The samples created or tested with SOUPE method can be plotted on the
phylogeny of the master-group using the package phytools (Revell, 2012). Here is
a code to reformat the output data of SOUPE functions to use phytools functions. 


``` r
tree.species <- ex_tree$tip.label

#creates a dataframe to classify if the species are present in the sample or not
presabs <- data.frame((matrix(NA,ncol = 2, nrow = length(tree.species)))) 
colnames(presabs) <- c("Species", "Presence")

presence <- tree.species %in% ex_tested_sample

for (i in 1:length(tree.species)){
  presabs$Species[[i]] <- tree.species[i]
  if (presence[[i]] == TRUE) {
    presabs$Presence[[i]] <- c("present")
    }
  else {
    presabs$Presence[[i]] <- c("absent")
  }
}

color = c("cyan3", "white")
trait <- presabs$Presence
names(trait) <- presabs$Species
colors <- setNames(color, c("present","absent"))


#rectangular tree using phytools::dotTree
phytools::dotTree(tree,
                  trait,
                  colors = colors,
                  ftype = "i",
                  legend = FALSE)
                  
legend(x = "bottomleft",
       legend = c("present", "absent"),
       fill = colors,
       cex = 1,
       bty = "n",
       horiz = TRUE)

#fan tree using phytools::plotFanTree.wTraits
colors.list <- list(setNames(color, c("present","absent")))

phytools::plotFanTree.wTraits(tree,
                              trait,
                              lwd = 12,
                              colors = colors.list,
                              ftype = "off")
                              
legend(x = "bottomleft",
       legend = c("present", "absent"),
       fill = colors,
       cex = 1,
       bty = "n",
       horiz = TRUE)

```

**Advice:**

**- The sample might be plotted on a phylogeny that doesn't include all the species
  of the master-group (e.g., phylogeny of available species). In this case, keep
  in mind that the visualized distribution can be biased, especially in case
  where the sub-phylogeny is not representative of the master-group's phylogeny.
  The distribution of the sample on the tree may not represent its distribution
  relatively to the master-group. To ensure the representativeness of a sample
  relatively to its master-group, functions of the a posteriori process should
  be used.**
  


## Bibliography

Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity.
Biological Conservation, 61(1), 1–10.
https://doi.org/10.1016/0006-3207(92)91201-3

Revell, L. J. (2012). Phytools: An r package for phylogenetic comparative
biology (And other things). Methods in Ecology and Evolution, 3(2), 217‑223.
https://doi.org/10.1111/j.2041-210X.2011.00169.x

