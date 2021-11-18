## ---------------------------------------------------------------------------------------------------------------
##
## Extracting and clustering colours
## 
## File:    extractColour.R
## 
## Description:
## 
##    Imports an image and extracts pixel level colours to be used
##    in various applications of k-means clustering
## 
## ---------------------------------------------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------------------
## Packages ####

library(cluster)
library(factoextra)
library(data.table)
library(doParallel)
library(dplyr)
library(foreach)
library(grid)
library(imager)
library(magick)
library(pacman)
library(stringr)

options(scipen = 999)

## ---------------------------------------------------------------------------------------------------------------
## Parameters ####

#Image
imgName <- "Norway"
imgType <- "jpg"
imgNameFull <- paste0(imgName, ".", imgType)

#Output
outName <- imgName
outType <- "jpg"
outLoc <- paste0(imgName, "/", "R", "/")

#Create output locations
dir.create(imgName, showWarnings = F) #Doesn't crash if the directory already exists
dir.create(outLoc, showWarnings = F) #Doesn't crash if the directory already exists

#K-means
maxK <- 30
iterK <- 3
maxIter <- 31

#Parellel cluster
# cores <- detectCores()
# cl <- makeCluster(cores-1)
# registerDoParallel(cl)


## ---------------------------------------------------------------------------------------------------------------
## Open image and get pixel colour data

img <- dcast(setDT(as.data.frame(load.image(imgNameFull))), x + y ~ cc)
colnames(img) <- c("x", "y", "R", "G", "B")

coords <- img[, .(x, y)]
colour <- img[, .(R, G, B)]

nrow(unique(colour))

#We need to reduce the size of the data for hierarchical clustering to work
scale <- 0.05

colour2 <- colour %>%
  mutate_all(function(x) { round(x / scale) * scale}) %>%
  unique(.)

#Distance matrix
d <- dist(colour2, method = "euclidean")

#perform hierarchical clustering using Ward's method
final_clust_ward <- hclust(d, method = "ward.D2")
final_clust_avg  <- hclust(d, method = "average")
final_clust_sing <- hclust(d, method = "single")
final_clust_comp <- hclust(d, method = "complete")

#We will use the above output cutting the tree for differing number of clusters

head(coords)
head(colour)

centres <- data.table()
clusters <- data.table()


## ---------------------------------------------------------------------------------------------------------------
## anim
##
## Description:
## 
##    Takes an existing data.frame of pixel colour data (three columns
##    representing fractional R, G & B values) and applies kmeans for an 
##    increasing value of k up to the value passed
## ---------------------------------------------------------------------------------------------------------------
anim <- function(maxk) {
  for (k in seq(1:maxk)) {
    #Run k-means on the colour data
    #Had issues with the defaults so using the
    #MacQueen algorithm with up to 100 iterations
    model <- kmeans(colour, k, iter.max=100, algorithm="MacQueen")
    saveRDS(model, paste0(outLoc, "model_", k, ".RDS"))
    
    #Also do a copy on the rounded colour values
    #Possibly a better comparison to the hierarchical version
    model2 <- kmeans(colour2, k, iter.max=100, algorithm="MacQueen")
    saveRDS(model2, paste0(outLoc, "model2_", k, ".RDS"))

    #Take the model output and assign each co-ordinate
    #the colour of it's cluster. Then reform into a
    #dataframe that can be converted to the cimg class
    newimg <- setDT(cbind(coords, model$centers[model$cluster, ])) %>%
      melt(id.vars = c("x", "y")) %>%
      mutate(cc = case_when( #colour channels
        variable == "R" ~ 1,
        variable == "G" ~ 2,
        variable == "B" ~ 3
      )) %>%
      select(x, y, cc, value) %>%
      as.cimg(dim=c(max(coords$x), max(coords$y), 1, 3))

    #Save a copy of the k-means cluster output
    save.image(newimg, paste0(outLoc, outName, "_anim", k, "_R_kmeans.", outType), quality=1)
    
    #We need to stick the reduced model back onto the original image
    clust <- cbind(colour2, cluster = model2$cluster)
    
    clust2 <- clust %>%
      group_by(cluster) %>%
      summarise_all(mean)
    
    #Merge back to original
    colour3 <- colour %>%
      mutate(R2 = round(R/0.05)*0.05,
             G2 = round(G/0.05)*0.05,
             B2 = round(B/0.05)*0.05) %>%
      left_join(clust, by=c("R2"="R", "G2"="G", "B2"="B")) %>%
      select(-c("R", "R2","G","G2","B","B2")) %>%
      left_join(clust2,by="cluster") %>%
      select(-cluster)
    
    newimg2 <- setDT(cbind(coords, colour3)) %>%
      melt(id.vars = c("x", "y")) %>%
      mutate(cc = case_when( #colour channels
        variable == "R" ~ 1,
        variable == "G" ~ 2,
        variable == "B" ~ 3
      )) %>%
      select(x, y, cc, value) %>%
      as.cimg(dim=c(max(coords$x), max(coords$y), 1, 3))
    
    #Save a copy of the k-means cluster output
    save.image(newimg2, paste0(outLoc, outName, "_anim", k, "_R_kmeans2.", outType), quality=1)

    #Next create an extension image to show the extracted colours
    #The background will be grey
    #This may cause issue if grey is an extracted colour
    extension <- expand.grid(x=seq(1:max(coords$x)),
                             y=seq(1:(max(coords$y)/4)),
                             #Grey background
                             R=0.6,
                             G=0.6,
                             B=0.6)
    
    extension2 <- extension

    #Keep track of the cluster centres labelling them as having k centres and numbering them 1 to k
    centres <- setDT(as.data.frame(model$centers))[, k:=k][, c:=seq(1:.N)][, method:="kmeans"]
    centres <- rbind(centres, setDT(as.data.frame(model2$centers))[, k:=k][, c:=seq(1:.N)][, method:="kmeans2"])
    
    #Also keep track of the cluster assigned to each pixel - probably allows us to re-create an image later if we want to
    clusters <- setDT(data.frame(cluster=model$cluster))[, k:=k][, method:="kmeans"]
    clusters <- rbind(clusters,setDT(data.frame(cluster=model2$cluster))[, k:=k][, method:="kmeans2"])
    
    #Some logic to add in boxes for each colour extracted
    xMargin <- 10
    yMargin <- 25
    width <- (max(coords$x) - ((ceiling(maxk/2) + 1) * xMargin))/ceiling(maxk/2)
    height <- (max(coords$y)/4 - (yMargin *3))/2
    for (c in seq(1:k)) {
      cMod = ifelse(c%%ceiling(maxk/2)==0, maxk/2, c%%ceiling(maxk/2))
      xMin <- (cMod - 1) * (xMargin + width) + xMargin + 1
      xMax <- cMod * (xMargin + width)
      yMin <- ifelse(c <= ceiling(maxk/2), yMargin, yMargin * 2 + height) + 1
      yMax <- ifelse(c <= ceiling(maxk/2), yMargin + height, yMargin * 2 + height * 2)

      extension <- extension %>%
        mutate(R = case_when(
                    x > xMin & x < xMax & y > yMin & y < yMax ~ model$centers[c, "R"],
                    TRUE ~ R),
               G = case_when(
                    x > xMin & x < xMax & y > yMin & y < yMax ~ model$centers[c, "G"],
                    TRUE ~ G),
               B = case_when(
                    x > xMin & x < xMax & y > yMin & y < yMax ~ model$centers[c, "B"],
                    TRUE ~ B)
        )
      
      extension2 <- extension2 %>%
        mutate(R = case_when(
                    x > xMin & x < xMax & y > yMin & y < yMax ~ model2$centers[c, "R"],
                    TRUE ~ R),
               G = case_when(
                    x > xMin & x < xMax & y > yMin & y < yMax ~ model2$centers[c, "G"],
                    TRUE ~ G),
               B = case_when(
                    x > xMin & x < xMax & y > yMin & y < yMax ~ model2$centers[c, "B"],
                    TRUE ~ B)
        )
    }

    extension <- setDT(extension)  %>%
      melt(id.vars = c("x", "y")) %>%
      mutate(cc = case_when( #colour channels
        variable == "R" ~ 1,
        variable == "G" ~ 2,
        variable == "B" ~ 3
      )) %>%
      select(x, y, cc, value) %>%
      as.cimg(dim=c(max(coords$x), max(coords$y)/4, 1, 3))
    
    extension2 <- setDT(extension2)  %>%
      melt(id.vars = c("x", "y")) %>%
      mutate(cc = case_when( #colour channels
        variable == "R" ~ 1,
        variable == "G" ~ 2,
        variable == "B" ~ 3
      )) %>%
      select(x, y, cc, value) %>%
      as.cimg(dim=c(max(coords$x), max(coords$y)/4, 1, 3))

    final <- list(newimg, extension) %>% imappend("y")
    final2 <- list(newimg2, extension2) %>% imappend("y")

    #save image with extensions
    save.image(final, paste0(outLoc, outName, "_anim", k, "_R_k_means_ext.", outType), quality=1)
    save.image(final2, paste0(outLoc, outName, "_anim", k, "_R_k_means2_ext.", outType), quality=1)
    
    #Hierarchical clusters
    clust <- cbind(colour2, cluster = cutree(final_clust_ward, k=k)) 
    
    #Get average colour for each cluster
    clust2 <- clust %>%
      group_by(cluster) %>%
      summarise_all(mean)
    
    #Merge back to original
    colour3 <- colour %>%
      mutate(R2 = round(R/0.05)*0.05,
             G2 = round(G/0.05)*0.05,
             B2 = round(B/0.05)*0.05) %>%
      left_join(clust, by=c("R2"="R", "G2"="G", "B2"="B")) %>%
      select(-c("R", "R2","G","G2","B","B2")) %>%
      left_join(clust2, by="cluster") 
    
    newimg <- setDT(cbind(coords, colour3)) %>%
      select(-cluster)%>%
      melt(id.vars = c("x", "y")) %>%
      mutate(cc = case_when( #colour channels
        variable == "R" ~ 1,
        variable == "G" ~ 2,
        variable == "B" ~ 3
      )) %>%
      select(x, y, cc, value) %>%
      as.cimg(dim=c(max(coords$x), max(coords$y), 1, 3))
    
    #Save a copy of the hierarchical cluster output
    save.image(newimg, paste0(outLoc, outName, "_anim", k, "_R_hier.", outType), quality=1)
    
    
    #Next create an extension image to show the extracted colours
    #The background will be grey
    #This may cause issue if grey is an extracted colour
    extension <- expand.grid(x=seq(1:max(coords$x)), 
                             y=seq(1:(max(coords$y)/4)), 
                             #Grey background
                             R=0.6, 
                             G=0.6, 
                             B=0.6)
    
    
    #Some logic to add in boxes for each colour extracted
    xMargin <- 10
    yMargin <- 25
    width <- (max(coords$x) - ((ceiling(maxk/2) + 1) * xMargin))/ceiling(maxk/2)
    height <- (max(coords$y)/4 - (yMargin *3))/2
    for (c in seq(1:k)) {
      cMod = ifelse(c%%ceiling(maxk/2)==0, maxk/2, c%%ceiling(maxk/2))
      xMin <- (cMod - 1) * (xMargin + width) + xMargin + 1
      xMax <- cMod * (xMargin + width)
      yMin <- ifelse(c <= ceiling(maxk/2), yMargin, yMargin * 2 + height) + 1
      yMax <- ifelse(c <= ceiling(maxk/2), yMargin + height, yMargin * 2 + height * 2)
      
      r <- as.numeric(clust2[c, "R"])
      g <- as.numeric(clust2[c, "G"])
      b <- as.numeric(clust2[c, "B"])
      
      extension <- extension %>%
        mutate(R = case_when(
          x > xMin & x < xMax & y > yMin & y < yMax ~ r,
          TRUE ~ R),
          G = case_when(
            x > xMin & x < xMax & y > yMin & y < yMax ~ g,
            TRUE ~ G),
          B = case_when(
            x > xMin & x < xMax & y > yMin & y < yMax ~ b,
            TRUE ~ B)
        )
    }
    
    extension <- setDT(extension)  %>%
      melt(id.vars = c("x", "y")) %>%
      mutate(cc = case_when( #colour channels
        variable == "R" ~ 1,
        variable == "G" ~ 2,
        variable == "B" ~ 3
      )) %>%
      select(x, y, cc, value) %>%
      as.cimg(dim=c(max(coords$x), max(coords$y)/4, 1, 3))
    
    final <- list(newimg, extension) %>% imappend("y")
    
    #save image with extensions
    save.image(final, paste0(outLoc, outName, "_anim", k, "_R_hier_ext.", outType), quality=1)
    
    
    #Keep track of the cluster centres labelling them as having k centres and numbering them 1 to k
    centres <- rbind(centres, setDT(clust2 %>% select(-cluster))[, k:=k][, c:=seq(1:.N)][, method:="hier"])
    
    #Also keep track of the cluster assigned to each pixel - probably allows us to re-create an image later if we want to
    clusters <- rbind(clusters,setDT(colour3 %>% select(cluster))[, k:=k][, method:="hier"])
  }
  return(list(centres, clusters))
}

output <- anim(maxK)

#Create start and end images

base <- cbind(coords, colour) %>%
  melt(id.vars = c("x", "y")) %>%
  mutate(cc = case_when( #colour channels
    variable == "R" ~ 1,
    variable == "G" ~ 2,
    variable == "B" ~ 3
  )) %>%
  select(x, y, cc, value) %>%
  as.cimg(dim=c(max(coords$x), max(coords$y), 1, 3))

extension <- expand.grid(x=seq(1:max(coords$x)), y=seq(1:(max(coords$y)/4)), R=0.6, G=0.6, B=0.6)
extension <- setDT(extension)  %>%
  melt(id.vars = c("x", "y")) %>%
  mutate(cc = case_when( #colour channels
    variable == "R" ~ 1,
    variable == "G" ~ 2,
    variable == "B" ~ 3
  )) %>%
  select(x, y, cc, value) %>%
  as.cimg(dim=c(max(coords$x), max(coords$y)/4, 1, 3))

final <- list(base, extension) %>% imappend("y")
save.image(base, paste0(outLoc, outName, "_anim_0_R_hier.", outType), quality=1)
save.image(final, paste0(outLoc, outName, "_anim_0_R_hier_ext.", outType), quality=1)
save.image(base, paste0(outLoc, outName, "_anim_0_R_kmeans.", outType), quality=1)
save.image(final, paste0(outLoc, outName, "_anim_0_R_kmeans_ext.", outType), quality=1)
save.image(base, paste0(outLoc, outName, "_anim_0_R_kmeans2.", outType), quality=1)
save.image(final, paste0(outLoc, outName, "_anim_0_R_kmeans2_ext.", outType), quality=1)


# model <- readRDS(paste0(outLoc, "model_", maxK, ".RDS"))
# centres <- rbind(centres, setDT(as.data.frame(model$centers))[, k:=maxK][, c:=seq(1:.N)])
# clusters <- rbind(clusters, setDT(as.data.frame(model$cluster))[, k:=maxK])


clust <- cbind(colour2, cluster= cutree(final_clust, k=maxK)) 

#Get average colour for each cluster
clust2 <- clust %>%
  group_by(cluster) %>%
  summarise_all(mean)

extension <- expand.grid(x=seq(1:max(coords$x)), y=seq(1:(max(coords$y)/4)), R=0.6, G=0.6, B=0.6)
#Some logic to add in boxes for each colour extracted
xMargin <- 10
yMargin <- 25
width <- (max(coords$x) - ((ceiling(maxK/2) + 1) * xMargin))/ceiling(maxK/2)
height <- (max(coords$y)/4 - (yMargin *3))/2
for (c in seq(1:maxK)) {
  cMod = ifelse(c%%ceiling(maxK/2)==0, maxK/2, c%%ceiling(maxK/2))
  xMin <- (cMod - 1) * (xMargin + width) + xMargin + 1
  xMax <- cMod * (xMargin + width)
  yMin <- ifelse(c <= ceiling(maxK/2), yMargin, yMargin * 2 + height) + 1
  yMax <- ifelse(c <= ceiling(maxK/2), yMargin + height, yMargin * 2 + height * 2)
  
  r <- as.numeric(clust2[c, "R"])
  g <- as.numeric(clust2[c, "G"])
  b <- as.numeric(clust2[c, "B"])
  
  extension <- extension %>%
    mutate(R = case_when(
      x > xMin & x < xMax & y > yMin & y < yMax ~ r,
      TRUE ~ R),
      G = case_when(
        x > xMin & x < xMax & y > yMin & y < yMax ~ g,
        TRUE ~ G),
      B = case_when(
        x > xMin & x < xMax & y > yMin & y < yMax ~ b,
        TRUE ~ B)
    )
}

extension <- setDT(extension)  %>%
  melt(id.vars = c("x", "y")) %>%
  mutate(cc = case_when( #colour channels
    variable == "R" ~ 1,
    variable == "G" ~ 2,
    variable == "B" ~ 3
  )) %>%
  select(x, y, cc, value) %>%
  as.cimg(dim=c(max(coords$x), max(coords$y)/4, 1, 3))

final <- list(base, extension) %>% imappend("y")
save.image(base, paste0(outLoc, outName, "_final_R_hier.", outType), quality=1)
save.image(final, paste0(outLoc, outName, "_final_R_hier_ext.", outType), quality=1)

createGIF <- function(pattern, fps, out) {
  list.files(path=outLoc, pattern = pattern, full.names = TRUE) %>% 
    str_sort(numeric=TRUE) %>%
    image_read() %>% # reads each path file
    image_join() %>% # joins image
    image_animate(fps=fps) %>% # animates, can opt for number of loops
    image_write(paste0(outLoc, outName, out)) # write to current dir
}

createGIF(pattern="*hier.jpg$", fps=2, out="_hier.gif")
createGIF(pattern="*kmeans.jpg$", fps=2, out="_kmeans.gif")
createGIF(pattern="*kmeans2.jpg$", fps=2, out="_kmeans2.gif")



colnames(output[[2]]) <- c("c", "k")
counts <- output[[2]][, .N, by=c("c", "k")][, perc:=N/max(N)][output[[1]], on=c("c", "k")][, col:=rgb(R, G, B)]

imgage <- jpeg::readJPEG(imgNameFull)


ggplot(counts[k==20], aes(x=c, y=N, fill=col)) +
  annotation_custom(rasterGrob(imgage,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")),
                    -Inf, Inf, -Inf, Inf) +
  geom_col() +
  scale_fill_identity()

counts[, variance:=apply(.SD, 1,var), .SDcols=c("R", "G", "B")]

distance <- copy(counts)[, t:=1]
distance <- distance[k==20, .(c, k, R, G, B, t)][distance[k==20, .(c, k, R, G, B, t)], on="t", allow.cartesian=T]

distance[, dist:= sqrt((R - i.R)^2 + (G - i.G)^2 + (B - i.B)^2)]
distance <- distance[dist != 0, lapply(.SD, mean), .SDcols="dist", by=c("c")]

counts <- counts[k==20][distance[, .(c, dist)], on="c"]

picks <- counts %>%
  arrange(desc(variance)) %>%
  slice(1:10) %>%
  arrange(desc(dist)) %>%
  slice(1:5)


ggplot(picks, aes(x=factor(c), y=N, fill=col)) +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc")), 
                    -Inf, Inf, -Inf, Inf) +
  geom_col() +
  scale_fill_identity() +
  ylim(c(0, max(picks$N)*5))


picks2 <- counts %>%
  arrange(desc(dist)) %>%
  slice(1:10) %>%
  arrange(desc(variance)) %>%
  slice(1:5)


ggplot(picks2, aes(x=factor(c), y=N, fill=col)) +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc")), 
                    -Inf, Inf, -Inf, Inf) +
  geom_col() +
  scale_fill_identity() +
  ylim(c(0, max(picks2$N)*5))


picks3 <- counts %>%
  mutate(dist*perc*variance) %>%
  slice(1:5)


ggplot(picks3, aes(x=factor(c), y=N, fill=col)) +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc")), 
                    -Inf, Inf, -Inf, Inf) +
  geom_col() +
  scale_fill_identity() +
  ylim(c(0, max(picks3$N)*5))

picks4 <- counts %>%
  mutate(pick = case_when(
    variance == max(variance) ~ 1,
    dist == max(dist) ~ 1,
    R*variance == max(R*variance) ~ 1,
    G*variance == max(G*variance) ~ 1,
    B*variance == max(B*variance) ~ 1,
    TRUE ~ 0
  )) %>%
  filter(pick == 1)


ggplot(picks4, aes(x=factor(c), y=N, fill=col)) +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc")), 
                    -Inf, Inf, -Inf, Inf) +
  geom_col() +
  scale_fill_identity() +
  ylim(c(0, max(picks4$N)*5))


picks_mod <- copy(picks)[, R:=min(1, R*1.1), by="c"]
picks_mod[, G:=min(1, G*1.1), by="c"]
picks_mod[, B:=min(1, B*1.1), by="c"]
picks_mod[, col:=rgb(R,G, B)]

ggplot(picks, aes(x=factor(c), y=N, fill=col)) +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc")), 
                    -Inf, Inf, -Inf, Inf) +
  geom_col() +
  scale_fill_identity() +
  ylim(c(0, max(picks$N)*5))

ggplot(picks_mod, aes(x=factor(c), y=N, fill=col)) +
  annotation_custom(rasterGrob(imgage, 
                               width = unit(1,"npc"), 
                               height = unit(1,"npc")), 
                    -Inf, Inf, -Inf, Inf) +
  geom_col() +
  scale_fill_identity() +
  ylim(c(0, max(picks_mod$N)*5))





picks_hsl <- picks %>%
  rowwise() %>%
  mutate(Cmax = max(R, G, B),
         Cmin = min(R, G, B),
         D = Cmax-Cmin,
         L = (Cmax + Cmin)/2,
         S = D / (1 - abs(2*L - 1)),
         H = case_when(
           D == 0 ~ 0,
           Cmax == R ~ 60 * mod((G - B)/D, 6),
           Cmax == G ~ 60 * (((B - R)/D) + 2),
           Cmax == B ~ 60 * (((R - G)/ D) + 4)),
         #Lighten
         L_mod = min(1, L*1.5),
         L_mod = L,
         #Now desaturate
         S_mod = min(1, S*0.7),
         S_mod = S,
         #Now change hue
         H_mod = mod(H + 350, 360),
         # H_mod = H,
         #And back to R, G, B
         C = (1 - abs(2*L_mod - 1)) * S_mod,
         X = C * (1 - abs(mod(H_mod/60, 2) - 1)),
         m = L_mod - C/2,
         R_mod = case_when(
           H_mod < 60 ~ C + m,
           H_mod < 120 ~ X + m,
           H_mod < 180 ~ 0 + m,
           H_mod < 240 ~ 0 + m,
           H_mod < 300 ~ X + m,
           TRUE ~ C + m),
         G_mod = case_when(
           H_mod < 60 ~ X + m,
           H_mod < 120 ~ C + m,
           H_mod < 180 ~ C + m,
           H_mod < 240 ~ X + m,
           H_mod < 300 ~ 0 + m,
           TRUE ~ 0 + m),
         B_mod = case_when(
           H_mod < 60 ~ 0 + m,
           H_mod < 120 ~ 0 + m,
           H_mod < 180 ~ X + m,
           H_mod < 240 ~ C + m,
           H_mod < 300 ~ C + m,
           TRUE ~ X + m),
         col2 = rgb(R_mod, G_mod, B_mod)
         )

ggplot(picks_hsl, aes(x=factor(c), y=N, fill=col)) +
  geom_col() +
  scale_fill_identity()

ggplot(picks_hsl, aes(x=factor(c), y=N, fill=col2)) +
  geom_col() +
  scale_fill_identity()




#Hierarchical
#define linkage methods
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
# ac <- function(x) {
#   agnes(colour2, method = x)$ac
# }

colour2 <- unique(colour %>%
                    mutate_all(function(x) { round(x / 0.05) * 0.05}))

# sapply(m, ac)

#perform hierarchical clustering using Ward's minimum variance
clust <- agnes(colour2, method = "ward")

#produce dendrogram
pltree(clust, cex = 0.6, hang = -1, main = "Dendrogram") 

#calculate gap statistic for each number of clusters (up to 10 clusters)
gap_stat <- clusGap(colour2, FUN = hcut, nstart = 15, K.max = 10, B = 30)

#produce plot of clusters vs. gap statistic
fviz_gap_stat(gap_stat)

#compute distance matrix
d <- dist(colour2, method = "euclidean")

#perform hierarchical clustering using Ward's method
final_clust <- hclust(d, method = "ward.D2" )

#cut the dendrogram into 4 clusters
groups <- cutree(final_clust, k=20)

#hierarchical centres
c <- cbind(colour2, cluster=groups) 

c2 <- c %>%
  group_by(cluster) %>%
  summarise_all(mean)

#Merge back to original
colour3 <- colour %>%
  mutate(R2 = round(R/0.05)*0.05,
         G2 = round(G/0.05)*0.05,
         B2 = round(B/0.05)*0.05) %>%
  left_join(c, by=c("R2"="R", "G2"="G", "B2"="B")) %>%
  select(-c("R", "R2","G","G2","B","B2")) %>%
  left_join(c2,by="cluster") %>%
  select(-cluster)

newimg <- setDT(cbind(coords, colour3)) %>%
  melt(id.vars = c("x", "y")) %>%
  mutate(cc = case_when( #colour channels
    variable == "R" ~ 1,
    variable == "G" ~ 2,
    variable == "B" ~ 3
  )) %>%
  select(x, y, cc, value) %>%
  as.cimg(dim=c(max(coords$x), max(coords$y), 1, 3))

ggplot() +
  annotation_custom(rasterGrob(newimg,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")),
                    -Inf, Inf, -Inf, Inf) 
