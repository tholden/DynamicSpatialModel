function CreateDefineNumbersMod( MaxNumber )
    FileID = fopen( 'DefineNumbers.mod', 'w' );
    if FileID > 0
        fprintf( FileID, '@#ifndef Numbers\n    @#define Numbers = [ ' );
        for i = 0 : MaxNumber
            fprintf( FileID, '"%d"', i );
            if i < MaxNumber
                fprintf( FileID, ', ' );
            end
        end
        fprintf( FileID, ' ]\n@#endif\n' );
        fclose( FileID );
    else
        disp( 'Error opening file. CreateDefineNumbersMod did not create DefineNumbers.mod.' );
    end
end
