
nt = 6210
ns = 1440*400
n.eof = 10 # extract only first EOFs



# connect to SciDB
library(scidb)
library(scidbst)
scidbconnect(host=Sys.getenv("scidb.host"),
             port = Sys.getenv("scidb.port"),
             username = Sys.getenv("scidb.user"), 
             password = Sys.getenv("scidb.pw"),
             auth_type = "digest",
             protocol = "https")


trmm.ref = scidb("TRMM3B42_DAILY_PREC")

# Create an array that maps time index to columns in the data matrix after removing incomplete images
t.complete =  unpack(
  subset(
    aggregate(
      subset(trmm.ref, "band1 >= 0"),
      by="t",
      FUN="count(*)"),
    "count = 576000"))$t

t.complete =  redimension(t.complete,schema="<i:int64>[t=0:6209,6209,0]")
t.complete = scidbeval(t.complete,name="TRMM3B42_NO_NA")
nt.nonempty =  aggregate(t.complete,FUN="count(*)")[]$count



# Build the data matrix
trmm.subset.X.ref =  
  redimension(
    transform(
      merge(trmm.ref, t.complete),
    s="int64(x*400+y)"),
  schema=paste("<band1:double>[i=0:", nt-1,", 32, 0, s=0:575999,32,0]",sep=""))

# detrend the data
trmm.subset.X.ref =  transform(
  merge(trmm.subset.X.ref , aggregate(
    trmm.subset.X.ref, by="s", "avg(band1)")), 
  prec_norm = "(band1 - band1_avg)")$prec_norm


# Run SVD to compute EOFs
trmm.subset.X.svd.R =  
  scidb::subarray(
    gesvd(trmm.subset.X.ref, type="right"),
    limits=c(0,0,n.eof-1,ns-1))

trmm.subset.X.svd.EOF.map =  
  redimension(
    transform(trmm.subset.X.svd.R ,y="int64(s % 400)", x="int64(floor(s / 400))"),
    schema=paste("<v:double>[y=0:399,400,0, x=0:1439,1440, 0, i=0:", n.eof - 1, ",1,0]", sep=""))


# Run all previous operations and store result as a new array
scidbeval(trmm.subset.X.svd.EOF.map,name="TRMM3B42_EOF")




