#Usage bash_andreaconcatpsl.sh out.psl

sed -i '1,5d' *.psl

echo "psLayout version 3" >> ${1}
echo "" >> ${1}
echo "match     mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts  tStarts" >> ${1}
echo "          match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count" >> ${1}
echo "---------------------------------------------------------------------------------------------------------------------------------------------------------------" >> ${1}

find -name "job*.psl" -exec tail -n +6 {} \; >> ${1}
