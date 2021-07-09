# echo "cleanup DotsConnection DotsConnection-00-00-05 in /workfs2/bes/mg20220135/Boss706a/Reconstruction"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtDotsConnectiontempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtDotsConnectiontempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=DotsConnection -version=DotsConnection-00-00-05 -path=/workfs2/bes/mg20220135/Boss706a/Reconstruction  $* >${cmtDotsConnectiontempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=DotsConnection -version=DotsConnection-00-00-05 -path=/workfs2/bes/mg20220135/Boss706a/Reconstruction  $* >${cmtDotsConnectiontempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtDotsConnectiontempfile}
  unset cmtDotsConnectiontempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtDotsConnectiontempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtDotsConnectiontempfile}
unset cmtDotsConnectiontempfile
return $cmtcleanupstatus

