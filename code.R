library(MASS)
library(car)
library(lme4)
library(HTSFilter)
library(parallel)
library(tidyverse)
library(data.table)

x<-read.delim("matrice_counts.txt",row.names="Geneid")

head(x)
names(x)[1] <- "QFT1"
names(x)[2] <- "QFT2"
names(x)[3] <- "QFT3"
names(x)[4] <- "QFT4"
names(x)[5] <- "QVT1"
names(x)[6] <- "QVT2"
names(x)[7] <- "QVT3"
names(x)[8] <- "QVT4"
names(x)[9] <- "QFL1"
names(x)[10] <- "QFL2"
names(x)[11] <- "QFL3"
names(x)[12] <- "QFL4"
names(x)[13] <- "QVL1"
names(x)[14] <- "QVL2"
names(x)[15] <- "QVL3"
names(x)[16] <- "QVL4"

names(x)[17] <- "WFT1"
names(x)[18] <- "WFT2"
names(x)[19] <- "WFT3"
names(x)[20] <- "WFT4"
names(x)[21] <- "WVT1"
names(x)[22] <- "WVT2"
names(x)[23] <- "WVT3"
names(x)[24] <- "WVT4"
names(x)[25] <- "WFL1"
names(x)[26] <- "WFL2"
names(x)[27] <- "WFL3"
names(x)[28] <- "WFL4"
names(x)[29] <- "WVL1"
names(x)[30] <- "WVL2"
names(x)[31] <- "WVL3"
names(x)[32] <- "WVL4"


ant_wide1 <- HTSBasicFilter(x,
                               method= "cpm",
                               cutoff.type = , # number of samples used for filtering
                               cutoff = , # CPM value
                               normalization=c("TMM"))

head(ant_wide1)


str(ant_wide1)

ant_wide1.df <- as.data.frame(ant_wide1$filteredData)
tibble::rownames_to_column(ant_wide1.df,var="Geneid") %>% as.data.table -> ant_wide1.dt

ant_tall = melt(ant_wide1.dt, id.vars=c("Geneid"))

head(ant_tall)


colnames(ant_tall) <- c("Geneid", "Replicate", "Count")
head(ant_tall)
ant_tall$Count <- round(ant_tall$Count, digits=0)
head(ant_tall)

##### merge together
ant_samples = fread('C:/Users/nmasm/Downloads/ant_metadata.txt')
tall_ant = merge(ant_tall, ant_samples, by = "Replicate")
tall_ant <- as.data.table(tall_ant)


##### make everything but count a factor
tall_ant$Replicate <- as.factor(tall_ant$Replicate)
tall_ant$Caste <- as.factor(tall_ant$Caste)
tall_ant$Habitat <- as.factor(tall_ant$Habitat)
tall_ant$Condition <- as.factor(tall_ant$Condition)




library(MASS)
library(car)
library(lme4)


run_tall_ant <- function(tall_ant, GID){
  print(GID)
  df1 <- tall_ant[Geneid == GID]
  glmm.res <- glmer(Count ~ Caste + Habitat + Condition + Caste:Habitat + Caste:Condition + Habitat:Condition + Caste:Habitat:Condition + (1 | Replicate), family = "poisson",
                    data = df1)
  fix <- car::Anova(glmm.res)
  aic <- AIC(glmm.res)
  results <- data.table("Geneid"= GID, 
                        "CasteB" = coef(summary(glmm.res))[2,1],
                        "CasteP" = coef(summary(glmm.res))[2,5],
                        "HabitatB" = coef(summary(glmm.res))[3,1],
                        "HabitatP" = coef(summary(glmm.res))[3,5],
                        "ConditionB" = coef(summary(glmm.res))[4,1],
                        "ConditionP" = coef(summary(glmm.res))[4,5],
                        "CasteHabitatB" = coef(summary(glmm.res))[5,1],
                        "CasteHabitatP" = coef(summary(glmm.res))[5,5],
                        "CasteConditionB" = coef(summary(glmm.res))[6,1],
                        "CasteConditionP" = coef(summary(glmm.res))[6,5],
                        "HabitatConditionB" = coef(summary(glmm.res))[7,1],
                        "HabitatConditionP" = coef(summary(glmm.res))[7,5],
                        "CasteHabitatConditionB" = coef(summary(glmm.res))[8,1],
                        "CasteHabitatConditionP" = coef(summary(glmm.res))[8,5],
                        "AIC" = aic[1])
  return(results)
}


run_tall_ant.safely = safely(run_tall_ant, otherwise = NULL)


results_parallel_tall_ant <- mclapply(unique(tall_ant$Geneid), function(GID){
  return(run_tall_ant.safely(tall_ant, GID))
}, mc.cores=detectCores())

results_parallel_tall_ant1 <- as.data.table(bind_rows(lapply(results_parallel_tall_ant, "[[", 1)))
fwrite(results_parallel_tall_ant1, '.csv')


results_parallel_tall_ant1$CasteHabitatadjust <- p.adjust(results_parallel_tall_ant1$CasteHabitatP, method = "BH")
results_parallel_tall_ant1$CasteConditionadjust <- p.adjust(results_parallel_tall_ant1$CasteConditionP, method = "BH")
results_parallel_tall_ant1$HabitatConditionadjust <- p.adjust(results_parallel_tall_ant1$HabitatConditionP, method = "BH")
results_parallel_tall_ant1$CasteHabitatConditionadjust <- p.adjust(results_parallel_tall_ant1$CasteHabitatConditionP, method = "BH")

results_parallel_tall_ant1$Casteadjust <- p.adjust(results_parallel_tall_ant1$CasteP, method = "BH")
results_parallel_tall_ant1$Habitatadjust <- p.adjust(results_parallel_tall_ant1$HabitatP, method = "BH")
results_parallel_tall_ant1$Conditionadjust <- p.adjust(results_parallel_tall_ant1$ConditionP, method = "BH")


results_parallel_tall_ant1three <- results_parallel_tall_ant1[CasteHabitatConditionadjust <= 0.05]
results_parallel_tall_ant1CH <- results_parallel_tall_ant1[CasteHabitatConditionadjust >= 0.05 & CasteHabitatadjust <= 0.05 & CasteConditionadjust >= 0.05 & HabitatConditionadjust >= 0.05]
results_parallel_tall_ant1CC <- results_parallel_tall_ant1[CasteHabitatConditionadjust >= 0.05 & CasteHabitatadjust >= 0.05 & CasteConditionadjust <= 0.05 & HabitatConditionadjust >= 0.05]
results_parallel_tall_ant1int <- results_parallel_tall_ant1[CasteHabitatConditionadjust >= 0.05 & CasteHabitatadjust >= 0.05 & CasteConditionadjust >= 0.05 & HabitatConditionadjust >= 0.05]
results_parallel_tall_ant1C <- results_parallel_tall_ant1int[Casteadjust < 0.05 & Habitatadjust > 0.05 & Conditionadjust > 0.05]
results_parallel_tall_ant1H <- results_parallel_tall_ant1int[Habitatadjust < 0.05 & Casteadjust > 0.05 & Conditionadjust > 0.05]
results_parallel_tall_ant1E <- results_parallel_tall_ant1int[Conditionadjust < 0.05 & Casteadjust > 0.05 & Habitatadjust > 0.05]


Here's a cleaner, more coherent version of your code that produces a single PCA plot:

```r
library(FactoMineR)
library(ggplot2)
library(dplyr)
library(data.table)
library(HTSFilter)
library(car)

# Load data
x <- read.delim("C:/Users/nmasm/Downloads/matrice_counts.txt", row.names = "Geneid")
ant_samples <- fread("C:/Users/nmasm/Downloads/ant_metadata.txt")

# Rename columns
colnames(x) <- c("QFT1", "QFT2", "QFT3", "QFT4", "QVT1", "QVT2", "QVT3", "QVT4",
                 "QFL1", "QFL2", "QFL3", "QFL4", "QVL1", "QVL2", "QVL3", "QVL4",
                 "WFT1", "WFT2", "WFT3", "WFT4", "WVT1", "WVT2", "WVT3", "WVT4",
                 "WFL1", "WFL2", "WFL3", "WFL4", "WVL1", "WVL2", "WVL3", "WVL4")

# Filter and normalize data
lauren_wide1 <- HTSBasicFilter(x, method = "cpm", cutoff.type = 16, cutoff = 1, normalization = "TMM")
lauren_wide1.df <- as.data.frame(lauren_wide1$filteredData)

# Log transform and transpose
x1 <- log(lauren_wide1.df + 1)
x1t <- t(x1)

# Perform PCA
res.pca.x <- PCA(x1t, scale.unit = TRUE, ncp = 3, graph = FALSE)
PCAres <- as.data.frame(res.pca.x$ind$coord)
PCAres <- tibble::rownames_to_column(PCAres, var = "Replicate")
PCAres.m <- merge(PCAres, ant_samples, by = "Replicate")
PCAres.m <- as.data.table(PCAres.m)

# Convert categorical variables to factors
PCAres.m <- PCAres.m %>% mutate(across(c(Replicate, Caste, Habitat, Condition), as.factor))

# Compute centroids
centroidsCa <- aggregate(cbind(Dim.1, Dim.2) ~ Caste, PCAres.m, mean)

# Define colors and shapes
color_vector <- ifelse(grepl("Labo", PCAres.m$Habitat),
                       ifelse(PCAres.m$Caste == "Queen", "#FF0000", "#0000FF"),
                       ifelse(PCAres.m$Caste == "Queen", "#FFA07A", "#87CEFA"))

shape_vector <- ifelse(PCAres.m$Caste == "Queen", 24, 22)

# Create plot
pdf(".pdf", width = 10, height = 8)
par(cex.lab = 1.45, cex.axis = 1.45)

plot(PCAres.m$Dim.1, PCAres.m$Dim.2, bg = color_vector, pch = shape_vector, cex = 1.5, 
     xlab = "PC1 (44.6%)", ylab = "PC2 (26.26%)", ylim = c(-200, 200), xlim = c(-150, 150))
points(centroidsCa$Dim.1, centroidsCa$Dim.2, pch = 3, col = "black", cex = 5)

dataEllipse(PCAres.m$Dim.1[PCAres.m$Caste == "Queen"], PCAres.m$Dim.2[PCAres.m$Caste == "Queen"], 
            col = "#FF0000", levels = c(0.95), center.pch = FALSE, plot.points = FALSE, add = TRUE)
dataEllipse(PCAres.m$Dim.1[PCAres.m$Caste == "Worker"], PCAres.m$Dim.2[PCAres.m$Caste == "Worker"], 
            col = "#0000FF", levels = c(0.95), center.pch = FALSE, plot.points = FALSE, add = TRUE)

mtext("PCA Plot", side = 3, line = -2, adj = 0.05, cex = 2)

library(dplyr)
library(data.table)
library(clusterProfiler)
library(GO.db)


CC <- fread("GO_cond_cas.txt")
#there was an extra blank 3rd column i had to remove manually first
anno <- fread("C:/Users/nmasm/Downloads/annotation_complete.csv", header = TRUE)
# all annotated genes have ".<number>.p1", but DE genes do not --> remove it?
anno$gene <- gsub(".[0-9].p1","",anno$gene)

# map go ids to term names
goterms <- Term(GOTERM)
# filter to those in annotation
anno_go <- goterms[anno$go]
# convert to df
anno_go_df = tibble(go = names(anno_go), description = anno_go)
# join to annotation
anno <- unique(left_join(anno, anno_go_df,by="go"))

# set up objects for clusterprofiler
go2gene=anno[, c("go", "gene")]
go2name=anno[, c("go", "description")]


# run enrichment
cc_enrichment = enricher(gene=CC$Geneid,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         minGSSize = 1,
                         maxGSSize = 1000000,
                         qvalueCutoff = 0.2,
                         TERM2GENE=go2gene,
                         TERM2NAME=go2name)

cc1 <- as.data.table(cc_enrichment)


