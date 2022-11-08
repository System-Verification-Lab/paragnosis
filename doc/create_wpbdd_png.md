To create a png from the compiled *asia* WPBDD or other representation, do the following:

    > pg compile asia --dot

This creates `.pg.output/asia.eps`. Convert that to `pdf` using `imagemagick`:

    > sudo apt-get install imagemagick
    > cd .pg.output
    > convert asia.eps asia.png

You now have a the png: `.pg.output/asia.png`.

## trouble shooting

If conversion returns with an error, saying permission denied, then remove the following from `/etc/ImageMagick-6/policy.xml`, and try again:

    <!-- disable ghostscript format types -->
    <policy domain="coder" rights="none" pattern="PS" />
    <policy domain="coder" rights="none" pattern="PS2" />
    <policy domain="coder" rights="none" pattern="PS3" />
    <policy domain="coder" rights="none" pattern="EPS" />
    <policy domain="coder" rights="none" pattern="PDF" />
    <policy domain="coder" rights="none" pattern="XPS" />
