
# connect to SciDB
library(scidb)
scidbconnect(host="localhost",
			 port = 8083,
			 username = "scidb", 
			 password = "xxxx.xxxx.xxxx",
			 auth_type = "digest",
			 protocol = "https")



# which datasets are in the database?
scidbst.ls()

# create a proxy object referencing the whole TRMM dataset
trmm.ref = scidbst("TRMM3B42_DAILY_PREC")
trmm.ref


# extract study region (South America)
trmm.ref = subarray(trmm.ref,textent(as.POSIXlt("1998-01-01"),as.POSIXlt("2015-01-01")), between=F)
trmm.ref.subregion = crop(trmm.ref,extent(-100, -25, -48, 20), between=F)



# remove images with missing values
array.dim = as.numeric(scidb_coordinate_bounds(trmm.ref.subregion)$length)
pixels.count = array.dim[1] * array.dim[2] 
pixels.count # pixels per image if there are no missing values
#2173


trmm.noinvalid = subset(trmm.ref.subregion, "band1 >= 0")     # leave out precipitation values < 0
#trmm.pixelcount = aggregate.sp(trmm.noinvalid,FUN="count(*)") # count available pixels per image  -> BUG

trmm.pixelcount = aggregate(trmm.noinvalid@proxy, by="t", FUN="count(*)") # count available pixels per image
t.complete = unpack(subset(trmm.pixelcount, paste("count > ", pixels.count -1 , sep="")))
t.complete@name
temp = t.complete[]
temp = trmm.pixelcount[]
iquery(trmm.ref.subregion@proxy@name, return = T)

head(trmm.ref@proxy)


# Create array that maps time index to columns in the data matrix after removing incomplete images
t.complete =  unpack(subset(aggregate(subset(trmm.subset.ref, "band1 >= 0"),by="t",FUN="count(*)"),"count = 576000"))$t
t.complete = redimension(t.complete,schema="<i:int64>[t=0:6209,6210,0]")
system.time(scidbeval(t.complete,name="TRMM3B42_SUBSET_COMPLETE_IMAGES_T")) # make result persistent
t.complete = scidb("TRMM3B42_SUBSET_COMPLETE_IMAGES_T")
nt.nonempty = aggregate(t.complete,FUN="count(*)")[]$count






## PCA, space as variable, subset 






#trmm.subset.ref = subset(subarray(trmm.ref$band1,limits = c(0,0,0,399,1439,364)), "band1 >= 0")
trmm.subset.ref = subarray(trmm.ref$band1,limits = c(0,0,0,ny-1,nx-1,nt-1))


# Create array that maps time index to columns in the data matrix after removing incomplete images
t.complete =  unpack(subset(aggregate(subset(trmm.subset.ref, "band1 >= 0"),by="t",FUN="count(*)"),"count = 576000"))$t
t.complete = redimension(t.complete,schema="<i:int64>[t=0:6209,6210,0]")
system.time(scidbeval(t.complete,name="TRMM3B42_SUBSET_COMPLETE_IMAGES_T")) # make result persistent
t.complete = scidb("TRMM3B42_SUBSET_COMPLETE_IMAGES_T")
nt.nonempty = aggregate(t.complete,FUN="count(*)")[]$count






# Build the data matrix
trmm.subset.X.ref = redimension(cast(transform( merge(trmm.subset.ref,t.complete),s="int64(x*400+y)"), "<band1:double,i:int64,s:int64> [y=0:399,2048,0,x=0:1439,2048,0,t=0:6209,1,0]"),schema="<band1:double>[i=0:6209, 32, 0, s=0:575999,32,0]") # surprisingly, an additional cast is required before redimension
trmm.subset.X.ref = subarray(trmm.subset.X.ref,c(0,0,nt.nonempty-1,ns-1))
system.time(scidbeval(trmm.subset.X.ref,name="TRMM3B42_SUBSET_EOF_S_X")) # 463 secs for one year






# detrend columns (subtract mean)
trmm.subset.X.ref = scidb("TRMM3B42_SUBSET_EOF_S_X")

# center and scale, check whether stdev is 0 for any pixel before!
#aggregate(subset(aggregate(trmm.subset.X.ref, by="s", FUN="stdev(band1)"),"band1_stdev = 0"), FUN="count(*)")[]
#trmm.subset.X.norm.ref = transform(merge(trmm.subset.X.ref , aggregate(trmm.subset.X.ref, by="s", "avg(band1),stdev(band1)")), prec_norm = "(band1 - band1_avg)/band1_stdev")$prec_norm

# only center
trmm.subset.X.norm.ref = transform(merge(trmm.subset.X.ref , aggregate(trmm.subset.X.ref, by="s", "avg(band1)")), prec_norm = "(band1 - band1_avg)")$prec_norm


# Run SVD to compute SVDs
trmm.subset.X.svd.R = gesvd(trmm.subset.X.norm.ref, type="right")
system.time(scidbeval(trmm.subset.X.svd.R,name="TRMM3B42_SUBSET_EOF_S_SVD_R")) # 497.78 secs for one year, 65610.713 ~ 18 h for whole dataset



# R is still transposed, so EOFs should be columns of R
trmm.subset.X.svd.R = scidb("TRMM3B42_SUBSET_EOF_S_SVD_R")
show(trmm.subset.X.svd.R)

n.eof = 10 # extract only first EOFs
trmm.subset.X.svd.EOF.map = redimension(transform(subarray(trmm.subset.X.svd.R,c(0,0,n.eof-1,ns-1)),y="int64(s % 400)", x="int64(floor(s / 400))"),schema=paste("<v:double NOT NULL>[y=0:399,32,0, x=0:1439, 32, 0, i=0:", n.eof - 1, ",1,0]", sep=""))
system.time(scidbeval(trmm.subset.X.svd.EOF.map,name="TRMM3B42_SUBSET_EOF_S")) # couple of secs for one year, couple of secs for whole dataset
trmm.subset.X.svd.EOF.map = scidb("TRMM3B42_SUBSET_EOF_S")




# extraction of i-th EOF
library(sp)
i = c(1,10)
eof.i =  subarray(trmm.subset.X.svd.EOF.map,c(0,0,i[1]-1,399,1439,i[2]-1))[]

a = cbind(sapply(i[1]:i[2],function(x) {
  eof.i$v[which(eof.i$i == x-1)]
}))
colnames(a) <- paste0("EOF.", i[1]:i[2] ,sep="")
eof = as.data.frame(a)
eof$x = -180+0.25*eof.i$x[which(eof.i$i == i[1]-1)]+0.25/2
eof$y = 50-0.25*eof.i$y[which(eof.i$i == i[1]-1)] -0.25/2




coordinates(eof) <- ~x+y
gridded(eof) = TRUE
proj4string(eof) <- CRS("+proj=longlat +datum=WGS84")
spplot(eof[c("EOF.1","EOF.2","EOF.3","EOF.4")],scales=list(TRUE))

spplot(eof, scales=list(TRUE))








