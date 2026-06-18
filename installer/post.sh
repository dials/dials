
mkdir -p ${PREFIX}/dialsbin
cp ${PREFIX}/bin/{dials,xia2,dxtbx,cctbx}* ${PREFIX}/dialsbin
cat >${PREFIX}/dials <<EOF
export PATH=${PREFIX}/dialsbin:\${PATH}
EOF
