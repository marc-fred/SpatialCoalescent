writeDataOutputInFile <- function(theta_sigma,theta_Y_k,theta_Y_r,theta_rate, file){
  # Writes arguments value and genetic data in file
  #
  # Args:
  #
  # Returns:
  # A file written 
  dir.create(path = paste(getwd(),"/Simulations", sep=""), showWarnings = FALSE)
  
  con <- file("Simulations/myFile", open = "w")
  writeLines(c("theta_sigma : ", as.character(theta_sigma)), con=con)
  writeLines(c("theta_Y_k : ", as.character(theta_Y_k)), con=con)
  writeLines(c("theta_Y_r : ", as.character(theta_Y_r)), con=con)
  writeLines(c("theta_rate : ", as.character(theta_rate)), con=con)
  writeLines(c("", "GENETICS:", ""), con=con)
  write.table(genetics, file=con, sep = "\t", quote = FALSE, col.names = FALSE, append=TRUE)
  close(con)
}
