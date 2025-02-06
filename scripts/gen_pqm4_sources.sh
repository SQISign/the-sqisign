#!/bin/bash

# This script should be run in the root folder of the repository, and creates pqm4 files in "src/pqm4/sqisign{1,3,5}"

if [ -d "src/pqm4" ]; then
    echo Destination folder src/pqm4 already exists. Delete it before running this script. Aborting.
    exit 1
fi

for LEVEL in 1 3 5
do
    LVL=lvl${LEVEL}
    DST_PATH=src/pqm4/sqisign_${LVL}/ref

    PQCGENKAT_SIGN_PQM4_BINARY=build/apps/PQCgenKAT_sign_pqm4_${LVL}

    if [ ! -f ${PQCGENKAT_SIGN_PQM4_BINARY} ]; then
        echo ${PQCGENKAT_SIGN_PQM4_BINARY} not found. Build it before running this script, or change the build folder in the script. Aborting.
        exit 1
    fi

    mkdir -p ${DST_PATH}

    # Run API generation script
    ${PQCGENKAT_SIGN_PQM4_BINARY}

    CPPFLAGS="-DRADIX_32 -DSQISIGN_BUILD_TYPE_REF -DSQISIGN_GF_IMPL_REF -DSQISIGN_VARIANT=${LVL} -DTARGET_ARM -DTARGET_OS_OTHER -DNDEBUG -DDISABLE_NAMESPACING -DBIG_PUBLIC_KEY_TESTS"
    PQM4_NAME="crypto_sign_sqisign_${LVL}_ref"

    echo "elf/${PQM4_NAME}_%.elf: CPPFLAGS+=${CPPFLAGS}" > ${DST_PATH}/config.mk
    echo "obj/lib${PQM4_NAME}.a: CPPFLAGS+=${CPPFLAGS}" >> ${DST_PATH}/config.mk

    cp include/{sig,sqisign_namespace}.h ${DST_PATH}/

    cp src/sqisign.c ${DST_PATH}/

    cp src/common/generic/include/{tools,tutil}.h ${DST_PATH}/

    cp src/ec/ref/lvlx/{basis,ec_jac,ec,isog_chains,xeval,xisog}.c ${DST_PATH}/
    cp src/ec/ref/include/{ec,isog}.h ${DST_PATH}/

    cp src/gf/ref/lvlx/{fp,fp2}.c ${DST_PATH}/
    cp src/gf/ref/include/{fp,fp2}.h ${DST_PATH}/

    cp src/hd/ref/lvlx/{hd.c,theta_isogenies.c,theta_isogenies.h,theta_structure.c,theta_structure.h} ${DST_PATH}/
    cp src/hd/ref/include/hd.h ${DST_PATH}/

    cp src/mp/ref/generic/mp.c ${DST_PATH}/
    cp src/mp/ref/generic/include/mp.h ${DST_PATH}/

    cp src/precomp/ref/${LVL}/include/{e0_basis,ec_params,encoded_sizes,fp_constants,hd_splitting_transforms}.h ${DST_PATH}/
    cp src/precomp/ref/${LVL}/{e0_basis,ec_params,hd_splitting_transforms}.c ${DST_PATH}/

    cp src/verification/ref/lvlx/{common,encode_verification,verify}.c ${DST_PATH}/
    cp src/verification/ref/include/verification.h ${DST_PATH}/
done

cp src/gf/ref/lvl1/fp_p5248_32.c src/pqm4/sqisign_lvl1/ref/
cp src/gf/ref/lvl3/fp_p65376_32.c src/pqm4/sqisign_lvl3/ref/
cp src/gf/ref/lvl5/fp_p27500_32.c src/pqm4/sqisign_lvl5/ref/
