
function list_alloc () {
cat ${files} | awk '{
i=index($0,"!")
if (i==0) {s0=$0} else {s0=substr($0,1,i-1)}
i=index(s0,"deallocate")
if (i==0) {
i=index(s0,"allocate")
if (i!=0) {
s=substr(s0,i)
gsub(" ","",s)
split(s,a,"(")
if (a[2] != "") {
s=a[2]; l=length(s)
if (substr(s,l) == ",") s=substr(s,1,l-1)
if (a[3] != "") print s
}}}}' - | sort | awk 'BEGIN{s=""}{
if ($1!=s && $1!="") {print $1; s=$1}
}' - >../alloc.lst
}
function list_dealloc () {

cat ${files} | awk -v d="${dir}" '{
i=index($0,"!")
if (i==0) {s0=$0} else {s0=substr($0,1,i-1)}
i=index(s0,"deallocate")
if (i!=0) {
s=substr(s0,i)
gsub(" ","",s)
i1=index(s,"("); i2=index(s,")")
s1=substr(s,i1+1,i2-i1-1)
if (i2!=0) print s1
}}' - | sort | awk 'BEGIN{s=""}{
if ($1!=s && $1!="") {print $1; s=$1}
}' - >../dealloc.lst
}
function compare_lists () {
diff alloc.lst dealloc.lst | awk -v d="${dir}" '
BEGIN{n=0;s=""}
{if ($1=="<") {na=na+1; sa=sa " \"" $2 "\""}
if ($1==">") {nd=nd+1; sd=sd " \"" $2 "\""}}
END{if (na!=0) print "Warning: variable" sa " not deallocated in",d
if (nd!=0) print "Warning: variable" sd " deallocated, but not allocated in",d}' -
}

BASE_DIR=$PWD
DIR_GLOBE="${BASE_DIR}/globe"
DIR_ATMOS="${BASE_DIR}/jam"
DIR_SOIL="${BASE_DIR}/nosoil"
DIR_VEG="${BASE_DIR}/jedi"

echo $DIR_GLOBE

  for dir in "${DIR_GLOBE}" "${DIR_ATMOS}" "${DIR_SOIL}" "${DIR_VEG}" ; do
    cd ${dir}
    files=`ls -1 *.*90`
    echo "Checking allocation and deallocation in directory ${dir}"
    list_alloc;
    list_dealloc;
    cd $BASE_DIR
    compare_lists;
    rm -f alloc.lst dealloc.lst
  done
