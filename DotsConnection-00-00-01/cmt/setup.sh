# echo "setup DotsConnection DotsConnection-00-00-05 in /workfs2/bes/mg20220135/Boss706a/Reconstruction"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtDotsConnectiontempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtDotsConnectiontempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=DotsConnection -version=DotsConnection-00-00-05 -path=/workfs2/bes/mg20220135/Boss706a/Reconstruction  -no_cleanup $* >${cmtDotsConnectiontempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=DotsConnection -version=DotsConnection-00-00-05 -path=/workfs2/bes/mg20220135/Boss706a/Reconstruction  -no_cleanup $* >${cmtDotsConnectiontempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtDotsConnectiontempfile}
  unset cmtDotsConnectiontempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtDotsConnectiontempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtDotsConnectiontempfile}
unset cmtDotsConnectiontempfile
return $cmtsetupstatus

