#!/usr/bin/env bash
sbatch <<EOT
#!/bin/bash

#SBATCH --output /n03data/ellien/Euclid_ERO/logs/head_%x.%j.out 
#SBATCH --error  /n03data/ellien/Euclid_ERO/logs/head_%x.%j.err
#SBATCH --time=$1
### This script works for any number of nodes, Ray will find and manage all resources
#SBATCH --nodes=1

### Give all resources on each node to a single Ray task, ray can manage the resources internally
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

# Load modules or your own conda environment here
conda init bash
source /home/ellien/.bashrc
conda activate dawis

head_node=$(hostname)
head_node_ip=$(hostname --ip-address)
# if we detect a space character in the head node IP, we'll
# convert it to an ipv4 address. This step is optional.
if [[ "$head_node_ip" == *" "* ]]; then
IFS=' ' read -ra ADDR <<<"$head_node_ip"
if [[ ${#ADDR[0]} -gt 16 ]]; then
  head_node_ip=${ADDR[1]}
else
  head_node_ip=${ADDR[0]}
fi
fi
port=6379

echo "STARTING HEAD at $head_node with time=$1"
echo "Head node IP: $head_node_ip"
ray start --head
sleep infinity