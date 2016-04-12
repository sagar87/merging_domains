#!/bin/zsh

#BSUB -q mpi
#BSUB -W 48:00
#BSUB -n 16
#BSUB -a openmp
#BSUB -o /cbscratch/hvoehri/scop_cluster/extract_a3m
#BSUB -R "span[hosts=1]"
#BSUB -R np16
#BSUB -R haswell
#BSUB -R cbscratch
#BSUB -J extract_a3m
#BSUB -m hh

#source paths.sh
source /etc/profile
source ~/.bashrc

mkdir -p /local/${USER}
MYLOCAL=$(mktemp -d --tmpdir=/local/${USER})

src_input=/cbscratch/hvoehri/scop_cluster/uniprot20_2012_10_1it_70
input_basename=$(basename ${src_input})
cp ${src_input}* ${MYLOCAL}
input=${MYLOCAL}/${input_basename}

a3m_database_extract -i ${input}_ca3m -o uniprot20_2012_10_1it_70_a3m -d ${input}_sequence -q ${input}_header

cp ${MYLOCAL}/uniprot20_2012_10_1it_70_a3m.ff{data,index} /cbscratch/hvoehri/scop_cluster/

#rm -f ${pdb70_build_dir}/pdb70_cs219.ff{data,index}