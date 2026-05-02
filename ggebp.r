# 
# ggebp.r
#
# Matias J. Curzi
#
# v1.3: 03/13/2015
#
# Performs the singular value decomposition (SVD) of a GxE matrix and constructs a GGE-biplot using the package ggplot2.
#
# The input for this function is a data frame of multi-environment trial (MET) data, in which genotypes
# are presented as rows and environments are presented as columns in a data matrix. The first column
# must be labeled 'Entry' and contain the genotype names. If the argument 'data' is not specified,
# a file selection window will appear that allows to choose a .csv file.
# 
# It is important to notice that this function does not admit missing values. If the dataset is unbalanced, rows
# or columns with missing data must be either removed or balanced using an appropriate imputation method.
#
# Alternatively, a dataset containing two factors and a response variable can be processed by this function.
# In this case, a column called 'Entry' must contain the genotype names, a second column called 'Env' must contain
# the environment names, and a third column called 'Yield' must have the response variable (Note: all possible 
# Entry x Env combinations must have a value). To use this input format, the argument 'mat' should be set to FALSE.
#
# An optional column called 'Group' is admitted to specify group names for the entries. This column is used
# to identify genotype groups by color in the biplot. In order to use this feature, the argument 'groups' must
# be set to TRUE.
#

ggebp <- function (data = NULL,  # A data frame object with MET data. If nothing is specified, a window will open to select a .csv file.
         env.cent = TRUE,        # Environment-centered biplot. Change to FALSE to scale variables (for an environment-standardized biplot).
         comps = c(1,2),         # Principal components to plot (PC1 and PC2 by default).
         output = NULL,          # If a name is specified between quotes, four .csv files will be created in the working directory.
         mat = TRUE,             # TRUE expects a data frame with entries in rows and environments in columns. FALSE expects a df with 3 columns: Entry, Env and Yield.
         mosaic = FALSE,         # Creates a mosaic plot showing the sum of squares partition in a new graphics device.
         groups = FALSE,         # If there are genotype groups, change to TRUE to identify them by color (an additional 'Group' column is required in the data set).
         title = "GGE-Biplot",   # Change biplot title (specify new title between quotes).
         obs.labels = FALSE,     # Change to TRUE to show observation (Entry) labels.
         var.labels = TRUE,      # Change to FALSE to hide variable (Environment) labels.
         angle = FALSE,          # Change to TRUE to assign environment labels the same angle as their vectors.
         obs.color = "black",    # Change fill color of entry symbols. Use "multi" to automatically set one color for each entry (ignored if groups = TRUE).
         var.color = "red4",     # Change line color of environment vectors. Use "multi" to automatically set one color for each environment.
         line.width = 0.5,       # Environment vector line width.
         obsname.size = 2.5,     # Entry labels text size.
         varname.size = 2.5,     # Environment labels text size.
         obs.factor = 1,         # Adjusts the observation coordinates to improve the plot appearance.
         var.factor = 1,         # Adjusts the variable coordinates to improve the plot appearance.
         sameLimits = TRUE,      # Change to FALSE to autoadjust x-axis and y-axis limits independently.
         xflip = TRUE,           # Change to FALSE to prevent autoflip of biplot x axis.
         yflip = FALSE,          # Change to TRUE to allow autoflip of biplot y axis.
         varname.adjust = -1.5) {  # Changes offset between variable vectors and their labels (works only if angle = TRUE). 


    library(ggplot2)

    if(is.null(data)) { 
        data <- read.csv(file.choose())
    }

    if(!mat) { 
        # Creates a list of dFrames containing entries and yields, one for each environment.
        datalist <- split(data[,c("Entry","Yield")], data$Env)

        # Saves the first element of the list in a data frame called 'dFrame'.
        names(datalist[[1]])[2] <- names(datalist)[[1]]
        dFrame <- datalist[[1]]

        # Loop to merge all the elements of the list in dFrame, creating a matrix
        # with genotypes in rows and environtmens in colums.
        for (i in 2:length(datalist)) { 
            names(datalist[[i]])[2] <- names(datalist)[[i]]
            dFrame <- cbind(dFrame, datalist[[i]][2])
        }
    } else {
        if("Group" %in% colnames(data)) {
	    dFrame <- subset(data, select = -Group)
        } else {
            dFrame <- data
        }
    }

    # Uses the Entry column as rownames and converts the data frame into a matrix.
    rownames(dFrame) <- dFrame$Entry
    dFrame <- as.matrix(subset(dFrame, select = -Entry ))

    # Environment-centered biplot vs. Environment-standardized biplot
    if(env.cent) {
	  Y <- scale(dFrame, scale=FALSE) 
    } else {
	  Y <- scale(dFrame, scale=TRUE) # Scales the variables (residuals divided by column standard deviation).
    }

########################################
# This section performs the SVD and calculates the partition of the total sum of squares:
# Formulas obtained from Laffont et al. 2007 (Crop Sci. 47:990–996) and Laffont et al. 2013 (Crop Sci. 53:2332–2341)

    num.gen <- nrow(Y)
    num.env <- ncol(Y)
    genMeans <- rowMeans(Y, na.rm=TRUE)

    YG <-  genMeans %*% t(rep(1, ncol(Y))) 
    YGE <- Y - YG 
    TSS <- sum(Y^2)                           # Total sum of squares
    SSG <- num.env * sum(genMeans^2)          # Sum of squares due to genotype
    SSGE <- sum(YGE^2)                        # Sum of squares due to G*E  

    Ysvd <- svd(Y)

  # Partition SSG and SSGE along each axis
    SSGk <- colSums(crossprod(YG, Ysvd$u)^2)
    SSGEk <- colSums(crossprod(YGE, Ysvd$u)^2)


  # Partition of SSG and SSGE along each axis as a percentage of the total.
    SSGkPerc <- round(SSGk/SSG*100,2)           # Percentage of SSG in axis k
    SSGEkPerc <- round(SSGEk/SSGE*100,2)        # Percentage of SSGE in axis k


  # Contribution of G and G*E effects within each axis:
    CGk <- round(SSGk/Ysvd$d^2*100, 2)
    CGEk <- round(SSGEk/Ysvd$d^2*100, 2)

########################################
# This section calculates observation and variable coordinates,
# adds the appropriate labels and genotype groups if present.

  # Getting observations and variables coordinates: 
  # (Adapted from function 'ggbiplot.r' by Vincent Q. Vu, 2011)
    nobs.factor <- sqrt(num.gen - 1) # Factor to ajust the observation scale.
    obs <- as.data.frame(nobs.factor * Ysvd$u) 
    env <- as.data.frame((1/nobs.factor) * Ysvd$v%*%diag(Ysvd$d))

  # Checks whether the requested axes are available:
    if(min(comps) < 1 | max(comps) > min(num.gen, num.env) | length(comps) > 2) {
        stop("There is an error in the PCs parameter.")
    }
    

  # Observations (genotypes)
    obsCoords <- obs[, comps] 
    obsCoords <- obsCoords * obs.factor

  # Variables (environments)    
    varCoords <- env[, comps]
    varCoords <- varCoords * var.factor
  
  # Axis names:
    names(obsCoords) <- c("coordX", "coordY")
    names(varCoords) <- names(obsCoords)

  # Flip the axes so that genotype ordinate is positively correlated with genotype means:
  # (Adapted from 'ggb' funtion by Laffont et. al, 2013)
    if(xflip & cor(genMeans, obsCoords$coordX) < 0) {
        varCoords$coordX <-  - varCoords$coordX
        obsCoords$coordX <-  - obsCoords$coordX
    }
    if(yflip & cor(genMeans, obsCoords$coordY) < 0) {
        varCoords$coordY <-  - varCoords$coordY
        obsCoords$coordY <-  - obsCoords$coordY
    } 

  # Observation and variable names assigend to a column:
    obsCoords$obsname <- rownames(dFrame)
    varCoords$varname <- colnames(dFrame)
  
  # Creates a grouping variable if required, used to identify genotype groups by color:
    if(groups) {
        obsCoords$group <- data$Group[match(obsCoords$obsname, data$Entry)]
    }

  # Variables for text label placement and angle calculation:
  # (Adapted from function 'ggbiplot.r' by Vincent Q. Vu, 2011)
    varCoords$angle <- with(varCoords, (180/pi) * atan(coordY / coordX))
    varCoords$hjust = with(varCoords, (1 + varname.adjust * sign(coordX)) / 2)

  # Axis labels with their proportion of explained variance:
    axis.labs <- paste('Axis', comps, sep='')
    axis.labs <- paste(axis.labs, sprintf('(%0.1f%%)', 100 * Ysvd$d[comps]^2 / sum(Ysvd$d^2)))


########################################
# Construction of the biplot using 'ggplot2' functions:
   
  # Observations:
    g <- ggplot() + 
        xlab(axis.labs[1]) + ylab(axis.labs[2]) +
        geom_hline(yintercept = 0, colour = "gray40", linetype= "dashed") +
        geom_vline(xintercept = 0, colour = "gray40", linetype= "dashed") +
        ggtitle(title) + theme(plot.title = element_text(vjust= 1)) 

    if(groups) {
        g <- g + geom_point(data = obsCoords, aes(x = coordX, y = coordY, label = obsname, fill = group), colour="black", shape=21)
    } else if (obs.color=="multi"){
        g <- g + geom_point(data = obsCoords, aes(x = coordX, y = coordY, label = obsname, fill = obsname), colour="black", shape=21)
    }else {
        g <- g + geom_point(data = obsCoords, aes(x = coordX, y = coordY, label = obsname), fill = obs.color, colour="black", shape=21)
    }

    # Observation labels:
    if(obs.labels) {
        if(groups) {  
            g <- g + geom_text(data = obsCoords, aes(x = coordX, y = coordY, label = obsname), colour = "black", vjust= -1, size = obsname.size)
        } else if (obs.color=="multi"){ 
            g <- g + geom_text(data = obsCoords, aes(x = coordX, y = coordY, label = obsname), colour = "black", vjust= -1, size = obsname.size)
        } else { 
            g <- g + geom_text(data = obsCoords, aes(x = coordX, y = coordY, label = obsname), colour = "black", vjust= -1, size = obsname.size)
        }
    }

    # Vectors (variables):
    if(var.color=="multi") {
        g <- g + geom_segment(data = varCoords, aes(x = 0, y = 0, xend = coordX, yend = coordY, colour = varname), size = line.width) +
             geom_point(data = varCoords, aes(x = coordX, y = coordY, colour = varname), size = 1.5)
    } else {
        g <- g + geom_segment(data = varCoords, aes(x = 0, y = 0, xend = coordX, yend = coordY), colour = var.color, size = line.width) +
             geom_point(data = varCoords, aes(x = coordX, y = coordY), colour = var.color, size = 1.5)   
    }
 
    # Vector labels:
    if(var.labels) {
        if(var.color=="multi") {
            if(angle) {
                g <- g + geom_text(data = varCoords, aes(label = varname, x = coordX, y = coordY, angle = angle, hjust = hjust, colour=varname), size = varname.size) 
            } else {
                g <- g + geom_text(data = varCoords, aes(x = coordX, y = coordY, label = varname, colour = varname), size = varname.size, vjust = 1.5) 
            }
        } else {
            if(angle) {
                g <- g + geom_text(data = varCoords, aes(label = varname, x = coordX, y = coordY, angle = angle, hjust = hjust), 
                    colour = var.color, size = varname.size) 
            } else {
                g <- g + geom_text(data = varCoords, aes(x = coordX, y = coordY, label = varname), colour = var.color, size = varname.size, vjust = 1.5) 
            }
        }
    }

    g <- g + scale_fill_discrete(name="Group") + scale_colour_discrete(guide=FALSE)

    if(sameLimits) {
      llim <- min(varCoords$coordX, obsCoords$coordX, varCoords$coordY, obsCoords$coordY)
      ulim <- max(varCoords$coordX, obsCoords$coordX, varCoords$coordY, obsCoords$coordY)
      g <- g + xlim(c(llim, ulim)) + ylim(c(llim, ulim))
    }

    if(!groups) {g <- g + theme(legend.position = "none")}

########################################
# Creates four tables: the partition of the Total SS, the proportion of SSG and SSE in each axis, the coordinates of the
# observations and the coordinates of the variables. These tables are returned by the function as a list object. Ff the
# 'output' parameter is specified, the tables are also saved as .csv files.

        PCs <- as.character()
        numrows <- min(nrow(Y), ncol(Y))

        PCs <- data.frame(PC = paste("PC", 1:numrows, sep=""))

        Percentage <- c("100 %", paste(round(SSG/TSS*100, 1), "%"), paste(round(SSGE/TSS*100, 1), "%")) 
        SSTable <- data.frame (SS = c("Total", "SSG", "SSGxE"), Value = c(TSS, SSG, SSGE), Percentage)    

        SSAxes <- data.frame(cbind(CGk, CGEk, SSGkPerc, SSGEkPerc))
        SSAxes <- cbind(PCs[1:numrows,], SSAxes[1:numrows,])
        colnames(SSAxes) <- c("PC", "%Contr_SSG_to_Axis",  "%Contr_SSGE_to_Axis", "%SSG", "%SSGE") 
        rownames(SSAxes) <- NULL
        colnames(obsCoords) <- c("x", "y", "Obs_Name")
        rownames(obsCoords) <- NULL
        colnames(varCoords)[1:3] <- c("x", "y", "Env_Name")
	rownames(varCoords) <- NULL

    if(!is.null(output)) {       
        write.csv(SSTable, file = paste(output, "_SS.csv", sep=""), row.names = FALSE)
        write.csv(SSAxes, file = paste(output, "_Axes.csv", sep=""), row.names = FALSE)
        write.csv(obsCoords, file = paste(output, "_ObsCoords.csv", sep=""), row.names = FALSE)
        write.csv(varCoords[,1:3], file = paste(output, "_VarCoords.csv", sep=""), row.names = FALSE)
    }

    result <- list(ssGGE = SSTable, ssAxes = SSAxes, obsCoords = obsCoords, varCoords = varCoords[,1:3])
    cat("\n")
    print(result$ssGGE)
    cat("\n")
    print(result$ssAxes)
    print(g) 


# Mosaic plot of SSG and SSGE partition:
    if(mosaic) {
      dev.new()
      mosdat <- data.frame(SSG =  SSGk, SSGE = SSGEk)
      mosdat <- as.matrix(mosdat)
      attr(mosdat, "dup.row.names") <- TRUE                            # Allows duplicate row names (also allows empty row names).
      rownames(mosdat) <- c("Axis1","Axis2", rep("", nrow(mosdat)-2))  # Only the first two axes have row names.
      mosaicplot(mosdat, main="", col=c("chartreuse4","orange2"), off=1)
      mtext("Mosaic plot of the TSS partition", line=1, cex=1.5)
    }


    cat("\nggebp: done!\n\n")
    return(invisible(result))  
}

########################################
#
# Changelog
# --------- 
#
# 01/26/2015: v1.0      Initial release in Spanish.
# 02/04/2015: v1.1      Solved file output error. Added SS partition in the screen output.
# 02/10/2015: v1.2      A list with useful data was added to the function output. Transtaled into English.
# 02/27/2015: v1.21     Added parameters to adjust cartesian coordinates and plot shape.
# 03/13/2015: v1.3      Use of 'svd' instead of 'prcomp'. Added mosaic plot of TSS. Added SS partition to console output.

