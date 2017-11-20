for /R %%f in (*.emf) do ( magick convert "%%f" "%%~npf.png" )
