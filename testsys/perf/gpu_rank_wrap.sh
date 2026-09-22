#!/bin/sh
# One GPU per MPI rank. rank k -> CUDA device k, by the LOCAL rank the launcher
# exports, so two ranks can never land on the same device and a 4-rank run is
# four A100s rather than four contexts on one.
#
# NO FALLBACK: if no launcher variable is present this is not an MPI rank, and
# defaulting to device 0 would put every rank on one GPU while the log said
# four. Refuse instead.
if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
  LR="$OMPI_COMM_WORLD_LOCAL_RANK"
elif [ -n "$MV2_COMM_WORLD_LOCAL_RANK" ]; then
  LR="$MV2_COMM_WORLD_LOCAL_RANK"
elif [ -n "$PMI_LOCAL_RANK" ]; then
  LR="$PMI_LOCAL_RANK"
else
  echo "gpu_rank_wrap.sh: no launcher local-rank variable in the environment;" >&2
  echo "cannot assign one GPU per rank. Refusing rather than defaulting to" >&2
  echo "device 0, which would run every rank on one device under an N-device label." >&2
  exit 3
fi
NGPU=$(nvidia-smi --query-gpu=index --format=csv,noheader | wc -l)
if [ "$LR" -ge "$NGPU" ]; then
  echo "gpu_rank_wrap.sh: local rank $LR but only $NGPU CUDA devices -- two ranks would share a device." >&2
  exit 3
fi
CUDA_VISIBLE_DEVICES="$LR"
export CUDA_VISIBLE_DEVICES
exec "$@"
