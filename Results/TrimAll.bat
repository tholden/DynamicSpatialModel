for /R %%f in (*.png) do ( magick convert "%%f" -trim "%%~npf-trim.png" )
