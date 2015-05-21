#!/bin/bash
cur_dir=`dirname $0`
pkg_dir=${cur_dir}/../pkg/
bin_dir=${cur_dir}/../bin/
mkdir -p ${pkg_dir}
mkdir -p ${bin_dir}

echo "downloading Star..."
star_url=https://github.com/alexdobin/STAR/archive/STAR_2.4.1d.tar.gz
star_tarball=${pkg_dir}${star_url##*/}
wget -P ${pkg_dir} --no-check-certificate -N ${star_url}
tar -zxvf ${star_tarball} -C ${pkg_dir}
star_dir=${pkg_dir}STAR-${star_url##*/}
star_dir=${star_dir%.tar.gz}
echo "copying star executable to bin/"
cp ${star_dir}/bin/Linux_x86_64/* ${bin_dir}
echo "export PATH..."
export LD_LIBRARY_PATH=${lib_dir}:$LD_LIBRARY_PATH
export PATH=${bin_dir}:$PATH
echo "done"
