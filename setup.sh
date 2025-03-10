source /nfs/data/1/yuhw/setup-sbnd.sh
source /nfs/data/1/yuhw/wire-cell-python/venv/bin/activate

path-prepend /nfs/data/1/yuhw/opt/lib LD_LIBRARY_PATH
path-prepend /nfs/data/1/yuhw/opt/bin PATH

WCT_SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
path-prepend ${WCT_SRC}/cfg/ WIRECELL_PATH
path-prepend ${WCT_SRC}/test/ BATS_LIB_PATH
path-prepend /nfs/data/1/yuhw/wire-cell-data WIRECELL_PATH

rs
export PS1="(app)"$PS1
