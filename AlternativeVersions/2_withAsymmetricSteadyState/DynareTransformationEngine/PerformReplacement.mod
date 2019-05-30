@#define OutputString = ""
@#define Continue = 1
@#for InputIndex in 1 : SpatialFunctionLength
    @#if Continue
        @#if InputString[InputIndex] == "#"
            @#define Continue = 0
        @#else
            @#if InputString[InputIndex] == "@"
                @#define OutputString = OutputString + ReplacementString
            @#else
                @#define OutputString = OutputString + InputString[InputIndex]
            @#endif
        @#endif
    @#endif
@#endfor
