# Alex on the 17/03/2020
# Function extraction

SavingInfo = "Extraction(x = my model)\n
Output: \n Vector of all the simple predictors in the model \n
Example : 
Extraction(MnbZT80.60.Gbis4)\n"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak))

Extraction.R <- function(x){
  print(paste0(getCall(x)[2]))
  A <- paste0(getCall(x)[2])
  A <- unlist(strsplit(A, "~",fixed = T))[2]
  A <- unlist(strsplit(A, "+",fixed = T))
  A <- grep(A,pattern = "|",fixed=T,value=T,invert=T)
  A <- unlist(strsplit(A, "*",fixed = T))
  A <- grep(A,pattern = ":",fixed=T,value=T,invert=T)
  A <- grep(A,pattern = "^[ ]$",value=T,invert=T) #Remove empty characters
  A <- sub("I(","", A, ignore.case = FALSE,fixed = T) # quadratique
  A <- sub("^2)","", A, ignore.case = FALSE,fixed = T) # quadratique
  A <- sub("offset(log(","", A, ignore.case = FALSE,fixed = T)
  A <- sub("))","", A, ignore.case = FALSE,fixed = T)
  A <- sub("\n ","", A, ignore.case = FALSE,fixed = T) # New line
  A <- sub(" ","", A, ignore.case = FALSE,fixed = T)
  A <- sub(" ","", A, ignore.case = FALSE,fixed = T) #Twice
  A <- sub(" ","", A, ignore.case = FALSE,fixed = T) #Twice
  A <- sub(" ","", A, ignore.case = FALSE,fixed = T) #Twice
  A <- unique(A)
  A
}
