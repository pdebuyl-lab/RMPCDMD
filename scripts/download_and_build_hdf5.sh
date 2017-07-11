#!/bin/bash
set -x

HDF5_NAME=hdf5-1.8.17
HDF5_TARBALL=${HDF5_NAME}.tar.gz
HDF5_SRC_URL=https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/${HDF5_NAME}/src/${HDF5_TARBALL}
HDF5_SRC_MD5SUM=7d572f8f3b798a628b8245af0391a0ca
BASEPATH="$(pwd)"
if [ -z "${TMPDIR}" ] ; then
   TMPDIR="${BASEPATH}"
fi

verify_sum() {
   if which md5sum ; then
      echo "${HDF5_SRC_MD5SUM}  ${HDF5_TARBALL}" | md5sum --quiet -c
   else
      LOCAL_SUM=$(openssl md5 -r ${HDF5_TARBALL})
      [ "$LOCAL_SUM" = "${HDF5_SRC_MD5SUM} *${HDF5_TARBALL}" ]
   fi
}

verify_sum
if [ "$?" != "0" ] ; then
    if [ $(which wget) ]
    then
	wget "${HDF5_SRC_URL}"
    else
	curl -L -O "${HDF5_SRC_URL}"
    fi
fi

verify_sum && (cd "${TMPDIR}" ; tar zxf "${BASEPATH}/${HDF5_TARBALL}")
if [ "$?" = "0" ] ; then
    cd "${TMPDIR}/${HDF5_NAME}"
    ./configure --enable-fortran --enable-fortran2003 --disable-shared --prefix="${BASEPATH}/_${HDF5_NAME}"
    make ${MAKEFLAGS} && make install
fi
