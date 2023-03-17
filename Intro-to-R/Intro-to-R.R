#Intro to R Lesson
#14.03.2023

# My first command via script!
x <-  3
y <-  5
number <- x+y

#create a numeric vector and name it as variable glength, numeric values in Mbp
glength <- c(4.6, 3000, 50000)

#create a character vector, storing the corresponding species
species <- c('ecoli', 'human', 'corn')
species

#combine those vectors
combined <- c(glength,species)

#create character vector and store as expression for different animals
expression <- c('low', 'high', 'medium', 'high', 'low', 'medium', 'high')

#turn expression vector into a factor
expression <- factor(expression)

#EXERCISE: create vector and turn it into a factor
samplegroup <- c('CTL', 'CTL', 'CTL', 'KO', 'KO', 'KO', 'OE', 'OE', 'OE')
samplegroup <- factor(samplegroup)

#create a df from vectors
df <- data.frame(species,glength)

#EXERCISE: create a dataframe from titles and pages
titles <- c('Catch-22', 'Pride and Prejudice', 'Nineteen Eighty Four')
pages <- c(453,432,328)
favorite_books = data.frame(titles,pages)

#getting to know lists
list1 <- list(species,df,number)

#EXERCISE: create a new list
list2 <- list(species, glength, number)

#EXERCISE: mean, sort
mean_glength <- mean(glength)
mean_glength

test <- c(1, NA, 2, 3, NA, 4)
mean_test <- mean(test, na.rm=T)
mean_test

sort(glength, decreasing = T)

#functions
square_it <- function(x){
  square <- x*x
  return(square)
}

#read files to R
options(stringsAsFactors = FALSE) #str werden nicht zwangsweise in einen faktor umgelesen, sondern verbleiben in einem Vektor
metadata <- read.csv(file='data/mouse_exp_design.txt')

#EXERCISE
proj_summary <- read.table(file='data/project-summary.txt',header=T,row.names=1)

#EXERCISE
class(glength)
class(metadata)

summary(proj_summary)

length(samplegroup)

dim(proj_summary)

class(rownames(metadata))

length(colnames(proj_summary))

#PRACTICE EXERCISES I
#1
temp_conv <- function(temp_f){
  temp_c <-(temp_f-32)*5/9
  return(temp_c+273.15)
}
temp_conv(70)

#2
round(temp_conv(70),1)


#PART 2
#create a vector called age
age <- c(15, 22, 45, 52, 73, 81)
age [5] #only the 5th value
age[-5] #all values but the 5th

age[c(3,5,6)] #retrieve several values
#or an other way:
idx <- c(3,5,6) #create a vector of the elements of interest
age[idx]

age[1:4]#retrieve first 4 values ordered in the "right" direction
age[4:1]#retrieve first 4 values ordered in the "wrong" direction

#EXERCISE
alphabets <- c('C', 'D', 'X', 'L', 'F')
alphabets[c(1,2,5)]
alphabets[-3]
alphabets[5:1]

#Use of logical operators
age > 50
age > 50 | age < 18
age
age[age>50|age<18] #nested
idx <- age>50|age<18
idx
age[idx]

#Indexing using which
which(age>50|age<18)
age[which(age>50|age<18)]

idx_num <- which(age>50|age<18)
age[idx_num]

#Applied to factors
expression[expression=='high']

#EXERCISE
samplegroup[samplegroup!='KO']

#Re-leveling of factors
expression
str(expression)

expression <- factor(expression, levels = c('low', 'medium', 'high'))
str(expression)

#EXERCISE
samplegroup <- factor(samplegroup, levels = c('KO', 'CTL', 'OE'))
str(samplegroup)

#Install packages from CRAN (latest versions, repository for installation)
install.packages('ggplot2')

#Alternatively, packages can also be installed from Bioconductor
#install.packages('BiocManager') --> only needs to be installed once in R
#BiocManager::install("ggplot2")

#Also, installation from source is possible
#install.packages("~/Downloads/ggplot2_1.0.1.tar.gz, type='source', repos=NULL)

#Load the installed package into your environment, so it is ready for use
library(ggplot2)


#extract values from df
metadata[1,1]
metadata[1,3]
metadata[,3]

#preventing to drop from a df structure to the easiest possible structure in R, in this case a vector
metadata[,3, drop=F]

# Dataframe containing first two columns
metadata[ , 1:2]

# Data frame containing first, third and sixth rows
metadata[c(1,3,6), ] 

# Extract the celltype column for the first three samples
metadata[c("sample1", "sample2", "sample3") , "celltype"] 

# Check column names of metadata data frame
colnames(metadata)

# Check row names of metadata data frame
rownames(metadata)

# Extract the genotype column
metadata$genotype 

# Extract the first five values/elements of the genotype column
metadata$genotype[1:5]

#EXERCISE
metadata[c('sample2','sample8'),c('genotype','replicate')]

metadata$replicate[c(4,9)]

metadata[,'replicate', drop=F]

metadata$celltype == 'typeA'

logical_idx <- metadata$celltype == 'typeA'

metadata[logical_idx,]

which(metadata$replicate>1)

idx <- which(metadata$replicate>1)
metadata[idx,]

#nesting
metadata[which(metadata$replicate>1), ]

#EXERCISE
metadata[which(metadata$genotype=='KO'),]

#Lists: to access an item you need to use double brackets!
list1[[2]]

#check class of list object
class(list1[[2]])

#access item in list item
list1[[2]][1,2]

#Advisable: do not mess with df itself, but write content to a variable and then alter it
#Using single brackets calling a list, it will return a list
list1[1]

#EXERCISE
random <- list(metadata,age,list1,samplegroup,number)
random[[4]]

#you can name list components
names(list1) <- c('species', 'df', 'number')
names(list1)

#extracting items from named lists via name
list1$df

#EXERCISE
names(random) <- c('metadata','age','list1','samplegroup','number')
random$age

rpkm_data <- read.csv('data/counts.rpkm.txt')
head(rpkm_data)

ncol(rpkm_data)
nrow(metadata)

#The %in% operator. Syntax: vector1 %in% vector2. Will test every item in vector 1 to be in vector 2. Do not have to have the same size
A <- c(1,3,5,7,9,11)   # odd numbers
B <- c(2,4,6,8,10,12)  # even numbers
# test to see if each of the elements of A is in B	
A %in% B

#RUN no2
A <- c(1,3,5,7,9,11)   # odd numbers
B <- c(2,4,6,8,1,5)  # add some odd numbers in 
# test to see if each of the elements of A is in B
A %in% B


intersection <- A %in% B
intersection
A[intersection]

#any/all function: is any/all value(s) from A in B?
any(A %in%B)
all(A %in%B)

#EXERCISE
any(B %in% A)
B[B %in% A]

A <- c(10,20,30,40,50)
B <- c(50,40,30,20,10)  # same numbers but backwards 

# test to see if each element of A is in B
A %in% B

# test to see if each element of A is in the same position in B
A == B

# use all() to check if they are a perfect match
all(A == B)

#Use functions on genetic data now
x <- rownames(metadata)
y <- colnames(rpkm_data)

all(x %in% y)
all(rownames(metadata) %in% colnames(rpkm_data))

x==y
all(x==y)

#EXERCISE
important_genes <- c("ENSMUSG00000083700", "ENSMUSG00000080990", "ENSMUSG00000065619", "ENSMUSG00000047945", "ENSMUSG00000081010", "ENSMUSG00000030970")

all(important_genes %in% rownames(rpkm_data))
extracted_rows = rpkm_data[c(rownames(rpkm_data) %in% important_genes),]
extracted_rows

#Matching and reordering
teaching_team <- c("Jihe", "Mary", "Meeta", "Radhika", "Will", "Emma")

teaching_team[c(2,4)]
reorder_teach <- teaching_team[c(5,4,6,2,1,3)]

#EXERCISE
first <- c("A","B","C","D","E")
second <- c("B","D","E","A","C")

for (x in c(1:5)){
  y <- second[c(first[[x]] == second)]
  second <- append(second, y)
}
second <- second[(length(second)/2+1):(length(second))]


#using the match() function. Input are always vectors!
first <- c("A","B","C","D","E")
second <- c("B","D","E","A","C")

match(first, second)
reorder_idx <- match(first, second)
second_reordered <- second[reorder_idx]


first <- c("A","B","C","D","E")
second <- c("D","B","A")  # remove values
match(first, second)
second[match(first, second)]

#NOTE: For values that don’t match by default return an NA value. You can specify what values you would have it assigned using nomatch argument. Also, if there is more than one matching value found only the first is reported


# Check row names of the metadata
rownames(metadata)

# Check the column names of the counts data
colnames(rpkm_data)

genomic_idx <- match(rownames(metadata), colnames(rpkm_data))
genomic_idx

# Reorder the counts data frame to have the sample names in the same order as the metadata data frame
rpkm_ordered  <- rpkm_data[ , genomic_idx]

all(rownames(metadata) == colnames(rpkm_ordered))

#EXERCISE
subset_rpkm <-  rpkm_ordered[,-c(2,9)]
subset_metadata <- metadata[match(colnames(subset_rpkm), rownames(metadata)),]

#Plotting and visualization
mean(rpkm_ordered$sample1)

#apply() and map() are function families which can be used instead of loops

library(purrr)

samplemeans <- map_dbl(rpkm_ordered, mean)

# Named vectors have a name assigned to each element instead of just referring to them as indices ([1], [2] and so on)
samplemeans

# Check length of the vector before adding it to the data frame
length(samplemeans)

# Create a numeric vector with ages. Note that there are 12 elements here
age_in_days <- c(40, 32, 38, 35, 41, 32, 34, 26, 28, 28, 30, 32)  


# Add the new vector as the last column to the new_metadata dataframe
new_metadata <- data.frame(metadata, samplemeans, age_in_days) 



#EXERCISE
animals <- read.csv(file='data/animals.txt')
class(animals)
dim(animals)

animals[animals$speed==40,]
animals[animals$color=='Tan',]
animals[animals$speed>50,'color', drop=F]
animals[animals$color=='Grey','color'] <- 'Gray'

animals_list <- list(animals$speed,animals$color)
names(animals_list) <- c('speed','color')

ctrl_samples <- data.frame(row.names = c("sample3", "sample10", "sample8", "sample4", "sample15"), date = c("01/13/2018", "03/15/2018", "01/13/2018", "09/20/2018","03/15/2018"))

length(which(rownames(ctrl_samples) %in% rownames(proj_summary)))
proj_summary_ctrl <- proj_summary[rownames(proj_summary) %in% rownames(ctrl_samples),]
#rownames(ctrl_samples) %in% rownames(proj_summary_ctrl)
m <- match(rownames(proj_summary_ctrl), rownames(ctrl_samples))
proj_summary_ctrl <- cbind(proj_summary_ctrl, batch=ctrl_samples[m,])

proj_summary_noctl <- proj_summary[which(proj_summary$treatment != "control"),]
keep <- map_lgl(proj_summary_noctl, is.numeric)
proj_summary_noctl <- proj_summary_noctl[,keep]








