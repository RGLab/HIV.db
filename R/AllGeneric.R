# TODO: Add comment
# 
# Author: mike
###############################################################################



setGeneric("isRoot", 
		function(object)
			standardGeneric("isRoot"))


setGeneric("getFrame", 
		function(object)
			standardGeneric("getFrame"))


setGeneric("getName", 
		function(object)
			standardGeneric("getName"))

setGeneric("getHIVdb", 
		function(object)
			standardGeneric("getHIVdb"))

setGeneric("getFeature", 
		function(object,...)
			standardGeneric("getFeature"))


setGeneric("getEpitope", 
		function(object,...)
			standardGeneric("getEpitope"))


setGeneric("getAnnotationTable", 
		function(object)
			standardGeneric("getAnnotationTable"))


setGeneric("getAntibodyTable", 
		function(object)
			standardGeneric("getAntibodyTable"))

setGeneric("getHXB2Coordinates", 
		function(object)
			standardGeneric("getHXB2Coordinates"))

setGeneric("HivFeatureList", 
		function(object,HIV_db)
			standardGeneric("HivFeatureList"))

setGeneric("EpitopeList", 
		function(object,HIV_db)
			standardGeneric("EpitopeList"))

setGeneric("getChildren", 
		function(object,...)
			standardGeneric("getChildren"))

setGeneric("getParent", 
		function(object,...)
			standardGeneric("getParent"))

setGeneric("getDNA", 
		function(object,...)
			standardGeneric("getDNA"))

setGeneric("getAA", 
		function(object,...)
			standardGeneric("getAA"))

setGeneric("getGenome",
		function(object,...)
			standardGeneric("getGenome"))

setGeneric("getRef",
		function(object)
			standardGeneric("getRef"))

