for /R %%f in (*.emf) do ( magick convert "%%f" -trim "%%~npf-trim.png" )
