
setGeneric(
  name = "ToStream",
  def = function(object, con, description) { return(standardGeneric("ToStream"))})

setGeneric(
  name = "myWrite",
  def = function(object, file) { return(standardGeneric("myWrite"))})

setGeneric(
  name = "getValues",
  def = function(object) { return(standardGeneric("getValues"))})

setGeneric(
  name = "myPlot",
  def = function(object) { return(standardGeneric("myPlot"))})

setGeneric(
  name = "getParameters",
  def = function(object) { return(standardGeneric("getParameters"))})

setGeneric(
  name = "applyFunction",
  def = function(object, xval) { return(standardGeneric("applyFunction"))})

setGeneric(
  name = "applyModel",
  def = function(object) { return(standardGeneric("applyModel"))})

setGeneric(
  name = "computeDistanceMatrix",
  def = function(object, con, description) { return(standardGeneric("computeDistanceMatrix"))})
