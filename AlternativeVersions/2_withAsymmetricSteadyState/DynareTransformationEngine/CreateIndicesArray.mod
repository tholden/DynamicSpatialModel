@#for CurrentIndex in 1 : SpatialPointsPerDimension
    @#define Indices = Indices + [ CurrentIndex ]
    @#define Depth = Depth + 1
    @#if Depth == SpatialDimensions
        @#define Name = ""
        @#for Dimension in 1 : SpatialDimensions
            @#define Name = Name + "_" + Numbers[Indices[Dimension]+1]
            @#define IndicesArray = IndicesArray + [ Indices[Dimension] ]
        @#endfor
        @#define IndicesStringArray = IndicesStringArray + [ Name ]
    @#else
        @#include "CreateIndicesArray.mod"
    @#endif
    @#define Depth = Depth - 1
    @#define Indices = [ Indices[ 1:Depth ] ]
@#endfor
