#!/bin/bash

DIR=$(mktemp -d)
BFILE=$(mktemp)
PROJECT=RMPCDMD
VERSION=$(git describe --always)
OUT="${DIR}/${PROJECT}-${VERSION}.tar"

git archive --prefix "${PROJECT}/"  HEAD -o "${OUT}"
git submodule foreach "git archive --prefix \"${PROJECT}/\${name}/\" HEAD -o \"${DIR}/\${name}.tar\" ; tar Af \"${OUT}\" \"${DIR}/\${name}.tar\""
mkdir -p "${DIR}/${PROJECT}/build"
(cd "${DIR}" ; mkdir -p "${PROJECT}"/build ; tar cf "${BFILE}" "${PROJECT}" ; tar Af "${OUT}" "${BFILE}" ; rm -r "${PROJECT}" "${BFILE}")

mv "${OUT}" .

rm -r ${DIR}

