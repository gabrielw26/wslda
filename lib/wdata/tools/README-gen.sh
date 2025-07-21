BINDIR=../bin/

echo "" > README.md # make empty file

echo '# wdata-copy' >> README.md
echo '```' >> README.md
$BINDIR/wdata-copy -h >> README.md
echo '```' >> README.md

echo '# wdata-subset' >> README.md
echo '```' >> README.md
$BINDIR/wdata-subset -h >> README.md
echo '```' >> README.md

echo '# wdata-section' >> README.md
echo '```' >> README.md
$BINDIR/wdata-section -h >> README.md
echo '```' >> README.md

echo '# wdata-append' >> README.md
echo '```' >> README.md
$BINDIR/wdata-append -h >> README.md
echo '```' >> README.md

echo '# wdata-datadim-up' >> README.md
echo '```' >> README.md
$BINDIR/wdata-datadim-up -h >> README.md
echo '```' >> README.md

echo '# wdata-interpolate' >> README.md
echo '```' >> README.md
$BINDIR/wdata-interpolate -h >> README.md
echo '```' >> README.md

echo '# wdata-d2f' >> README.md
echo '```' >> README.md
$BINDIR/wdata-d2f -h >> README.md
echo '```' >> README.md

echo '# wdata-f2d' >> README.md
echo '```' >> README.md
$BINDIR/wdata-f2d -h >> README.md
echo '```' >> README.md
