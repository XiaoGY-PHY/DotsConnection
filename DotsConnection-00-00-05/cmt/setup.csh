# echo "setup DotsConnection DotsConnection-00-00-05 in /workfs2/bes/mg20220135/Boss706a/Reconstruction"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtDotsConnectiontempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtDotsConnectiontempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=DotsConnection -version=DotsConnection-00-00-05 -path=/workfs2/bes/mg20220135/Boss706a/Reconstruction  -no_cleanup $* >${cmtDotsConnectiontempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=DotsConnection -version=DotsConnection-00-00-05 -path=/workfs2/bes/mg20220135/Boss706a/Reconstruction  -no_cleanup $* >${cmtDotsConnectiontempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtDotsConnectiontempfile}
  unset cmtDotsConnectiontempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtDotsConnectiontempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtDotsConnectiontempfile}
unset cmtDotsConnectiontempfile
exit $cmtsetupstatus

