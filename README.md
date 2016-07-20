# LCIOUtils
This is a small package which provides several convenience functions for interfacing with SLCIO files. These are

```jl
getpositions(filename, collectionName="PandoraPFOCollection", eventNum=1)
mappositions(filename, startevent, endevent, collectionname="PandoraPFOCollection")
mapcollections(f, filename, startevent, endevent, collections...)
mapparticles(f, filename, startevent, endevent, collectionname="PandoraPFOCollection")
maptruthparticles(f, filename, startevent, endevent)
```

In addition, some other convenience functions specific to this project are provided:
```jl
append(filename, contents) #append contents to file with given filename
inertiatensor(positions::Matrix) #compute the inertia tensor of the weighted points represented by the 4xn matrix `positions`
```
