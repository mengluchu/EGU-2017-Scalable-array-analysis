

# connect to SciDB
library(scidb)
library(scidbst)
scidbconnect(host=Sys.getenv("scidb.host"),
             port = Sys.getenv("scidb.port"),
             username = Sys.getenv("scidb.user"), 
             password = Sys.getenv("scidb.pw"),
             auth_type = "digest",
             protocol = "https")



# which datasets are in the database?
scidbst.ls()

scidbst.ls(extent = TRUE)

# create a proxy object referencing the whole TRMM dataset and
# show some metadata
trmm.ref = scidbst("TRMM3B42_DAILY_PREC")
trmm.ref


trmm.ref = subset(trmm.ref, "band1 >= 0")

# extract study region (South America)
trmm.ref.subregion = crop(trmm.ref,extent(-100, -25, -50, 20))
trmm.ref.subregion = subarray(trmm.ref.subregion,textent(as.POSIXlt("2005-01-01"),as.POSIXlt("2006-01-01")))
trmm.ref.subregion

# summarize pixel time series 
trmm.ref.subregion.summary = aggregate.t(trmm.ref.subregion,FUN="avg(band1)")

# run all previous queries and store result as a new array
scidbsteval(trmm.ref.subregion.summary, "TEMP")

# download result
require(raster)
result = as(scidbst("TEMP"), "RasterBrick")
plot(result)













