#pip
packname=`python -m pip wheel . --use-feature=in-tree-build --no-deps -w . | grep -o 'hifast-.*whl'`
echo "package name: "$packname
